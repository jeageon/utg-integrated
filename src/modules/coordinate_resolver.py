from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Any, Optional

from ..config import (
    ENSEMBL_LOOKUP,
    EBI_COORDINATES_URL,
    NCBI_ESEARCH,
    NCBI_ESUMMARY,
    UNIPROT_ENTRY_URL,
    UNIPROT_IDMAP_RUN,
    UNIPROT_IDMAP_RESULTS,
    UNIPROT_IDMAP_STATUS,
)
from ..models.data_schemas import GenomicCoordinates
from ..utils.api_client import ApiClient
from ..utils.coord_utils import apply_flank
from ..utils.exceptions import NoMappingError, ToolError


@dataclass
class ResolverResult:
    coordinates: GenomicCoordinates
    warnings: list[str]


_ENSEMBL_GENE_ID_RE = re.compile(r"^ENS\w*G\d+$", re.IGNORECASE)


def _normalize_ensembl_gene_id(raw: Any) -> Optional[str]:
    if not isinstance(raw, str):
        return None
    cleaned = raw.strip()
    if not cleaned:
        return None
    cleaned = cleaned.split(".", 1)[0]
    return cleaned if _ENSEMBL_GENE_ID_RE.match(cleaned) else None


def _extract_gene_from_ref_value(ref: Any) -> Optional[str]:
    if not isinstance(ref, dict):
        return None

    direct = ref.get("id")
    normalized = _normalize_ensembl_gene_id(direct)
    if normalized:
        return normalized

    for prop in ref.get("properties", []) or []:
        if not isinstance(prop, dict):
            continue
        key = str(prop.get("key") or "").lower()
        if key == "geneid":
            normalized = _normalize_ensembl_gene_id(prop.get("value"))
            if normalized:
                return normalized
    return None


class CoordinateResolver:
    def __init__(self, api_client: ApiClient) -> None:
        self.api = api_client

    def resolve(
        self,
        uniprot_id: str,
        flank_bp: int,
        flank_mode: str = "genomic",
        assembly_preference: str = "auto",
        taxid_filter: Optional[int] = None,
    ) -> ResolverResult:
        del assembly_preference
        warnings: list[str] = []

        entry_for_routing = self._fetch_uniprot_entry(uniprot_id)
        use_ncbi_first = self._is_bacterial_entry(entry_for_routing)
        if use_ncbi_first:
            warnings.append("microbial mode: attempting NCBI-first resolution")

        if use_ncbi_first:
            ncbi_lookup, ncbi_warnings = self._resolve_ncbi_gene(uniprot_id, uniprot_entry=entry_for_routing)
            warnings.extend(ncbi_warnings)
            if ncbi_lookup:
                return self._build_ncbi_result(uniprot_id, ncbi_lookup, flank_bp, flank_mode, taxid_filter, warnings)

        ensembl_gene_id, fallback_warnings = self._resolve_ensembl_gene(uniprot_id)
        warnings.extend(fallback_warnings)
        if ensembl_gene_id:
            lookup = self._lookup_ensembl_gene(ensembl_gene_id)
            if lookup:
                return self._build_ensembl_result(uniprot_id, lookup, flank_bp, flank_mode, taxid_filter, warnings)

        if not use_ncbi_first:
            ncbi_lookup, ncbi_warnings = self._resolve_ncbi_gene(uniprot_id, uniprot_entry=entry_for_routing)
            warnings.extend(ncbi_warnings)
            if not ncbi_lookup:
                detail = "; ".join(warnings) if warnings else "no valid mapping found"
                raise NoMappingError(f"No mapping found for {uniprot_id}: {detail}")
            return self._build_ncbi_result(uniprot_id, ncbi_lookup, flank_bp, flank_mode, taxid_filter, warnings)

        detail = "; ".join(warnings) if warnings else "no valid mapping found"
        raise NoMappingError(f"No mapping found for {uniprot_id}: {detail}")

    def _is_bacterial_entry(self, entry: Optional[dict[str, Any]]) -> bool:
        if not isinstance(entry, dict):
            return False
        organism = entry.get("organism") or {}
        if not isinstance(organism, dict):
            return False
        lineages = self._collect_lineages(organism)
        microbial_markers = {"bacteria", "archaea", "viral", "viruses"}
        if any(marker in lineages for marker in microbial_markers):
            return True
        organism_name = self._coerce_str(organism.get("scientificName") or organism.get("taxon-scientific-name"))
        return "bacteria" in organism_name.lower() or "archaea" in organism_name.lower()

    def _collect_lineages(self, organism: dict[str, Any]) -> set[str]:
        values: list[Any] = []
        for key in ("lineage", "lineages"):
            value = organism.get(key)
            if value is None:
                continue
            if isinstance(value, list):
                values.extend(value)
            else:
                values.append(value)
        normalized: set[str] = set()
        for item in values:
            if isinstance(item, dict):
                for key in ("scientificName", "name", "value", "taxon"):
                    candidate = self._coerce_str(item.get(key))
                    if candidate:
                        normalized.add(candidate.lower())
                continue
            candidate = self._coerce_str(item)
            if candidate:
                normalized.add(candidate.lower())
        return normalized

    def _build_ensembl_result(
        self,
        uniprot_id: str,
        lookup: dict[str, Any],
        flank_bp: int,
        flank_mode: str,
        taxid_filter: Optional[int],
        warnings: list[str],
    ) -> ResolverResult:
        gene_start = lookup["gene_start_1based"]
        gene_end = lookup["gene_end_1based"]
        if gene_start is None or gene_end is None:
            raise NoMappingError(f"Invalid coordinate span for {lookup.get('ensembl_gene_id', uniprot_id)}")
        if gene_start > gene_end:
            warnings.append(
                f"coordinate span was reordered for {lookup.get('ensembl_gene_id', uniprot_id)}: {gene_start}>{gene_end}"
            )
            gene_start, gene_end = gene_end, gene_start
        if taxid_filter is not None and lookup.get("taxid") != taxid_filter:
            warnings.append(f"taxid mismatch: got {lookup.get('taxid')}, requested {taxid_filter}")

        strand = lookup["strand"]
        ext_start, ext_end = apply_flank(gene_start, gene_end, flank_bp, flank_mode, strand)
        return ResolverResult(
            coordinates=GenomicCoordinates(
                uniprot_id=uniprot_id,
                ensembl_gene_id=lookup["ensembl_gene_id"],
                coordinate_source="ensembl",
                species=lookup["species"],
                assembly_name=lookup["assembly_name"],
                seq_region_name=lookup["seq_region_name"],
                gene_start_1based=gene_start,
                gene_end_1based=gene_end,
                strand=strand,
                display_name=lookup.get("display_name"),
                taxid=lookup.get("taxid"),
                ext_start_1based=ext_start,
                ext_end_1based=ext_end,
            ),
            warnings=warnings,
        )

    def _build_ncbi_result(
        self,
        uniprot_id: str,
        ncbi_lookup: dict[str, Any],
        flank_bp: int,
        flank_mode: str,
        taxid_filter: Optional[int],
        warnings: list[str],
    ) -> ResolverResult:
        gene_start = ncbi_lookup["gene_start_1based"]
        gene_end = ncbi_lookup["gene_end_1based"]
        if gene_start is None or gene_end is None:
            raise NoMappingError(f"Invalid coordinate span for {uniprot_id}")
        if gene_start > gene_end:
            warnings.append(f"coordinate span was reordered for {uniprot_id}: {gene_start}>{gene_end}")
            gene_start, gene_end = gene_end, gene_start
        if taxid_filter is not None and ncbi_lookup.get("taxid") is not None and ncbi_lookup.get("taxid") != taxid_filter:
            warnings.append(f"taxid mismatch: got {ncbi_lookup.get('taxid')}, requested {taxid_filter}")

        strand = ncbi_lookup["strand"]
        ext_start, ext_end = apply_flank(gene_start, gene_end, flank_bp, flank_mode, strand)
        return ResolverResult(
            coordinates=GenomicCoordinates(
                uniprot_id=uniprot_id,
                ensembl_gene_id=ncbi_lookup["ensembl_gene_id"],
                coordinate_source="ncbi",
                ncbi_accession=ncbi_lookup["ncbi_accession"],
                species=ncbi_lookup["species"],
                assembly_name=ncbi_lookup["assembly_name"],
                seq_region_name=ncbi_lookup["seq_region_name"],
                gene_start_1based=gene_start,
                gene_end_1based=gene_end,
                strand=strand,
                display_name=ncbi_lookup.get("display_name"),
                taxid=ncbi_lookup.get("taxid"),
                ext_start_1based=ext_start,
                ext_end_1based=ext_end,
            ),
            warnings=warnings,
        )

    def _resolve_ensembl_gene(self, uniprot_id: str) -> tuple[Optional[str], list[str]]:
        warnings: list[str] = []
        try:
            response = self.api.get(
                EBI_COORDINATES_URL.format(accession=uniprot_id),
                headers={"Accept": "application/json"},
            )
            payload = response.json_obj or []
            candidates: list[dict[str, Any]] = []
            if isinstance(payload, dict):
                payload = [payload]
            if not isinstance(payload, list):
                raise ToolError("Invalid EBI response structure")
            for item in payload:
                gene_id = self._extract_ensembl_gene_id(item)
                coords = self._parse_coords(item)
                if gene_id:
                    candidate = {"gene_id": gene_id, **coords}
                    candidates.append(candidate)
                if candidates:
                    if len(candidates) > 1:
                        warnings.append("multiple EBI coordinate candidates found; selected first valid one")
                    return candidates[0]["gene_id"], warnings
        except ToolError as exc:
            warnings.append(f"EBI coordinates lookup failed for {uniprot_id}: {exc}")

        try:
            return self._fallback_uniprot_mapping(uniprot_id)
        except ToolError as exc:
            warnings.append(f"UniProt mapping lookup failed for {uniprot_id}: {exc}")
            return None, warnings

    def _extract_ensembl_gene_id(self, entry: dict[str, Any]) -> Optional[str]:
        direct = entry.get("ensemblGeneId") or entry.get("ensembl_gene_id")
        normalized = _normalize_ensembl_gene_id(direct)
        if normalized:
            return normalized
        for ref in entry.get("crossReferences", []) or []:
            db_name = str(ref.get("dbDisplayName", "")).lower()
            if "ensembl" not in db_name:
                continue
            normalized = _extract_gene_from_ref_value(ref)
            if normalized:
                return normalized
        return None

    def _parse_coords(self, entry: dict[str, Any]) -> dict[str, Any]:
        loc = entry.get("genomicLocation") or entry.get("genomic_location") or {}
        start = loc.get("start") or entry.get("genomicStart") or entry.get("start")
        end = loc.get("end") or entry.get("genomicEnd") or entry.get("end")
        chrom = (
            loc.get("chromosome")
            or loc.get("seqRegion")
            or loc.get("seqRegionName")
            or entry.get("seqRegion")
            or entry.get("seq_region_name")
        )
        return {
            "start": _to_int(start),
            "end": _to_int(end),
            "chrom": str(chrom) if chrom is not None else None,
            "taxid": _to_int(entry.get("organism", {}).get("taxid") or entry.get("taxId") or entry.get("taxid")),
            "strand": _to_int(entry.get("strand")) or 1,
        }

    def _fallback_uniprot_mapping(self, uniprot_id: str) -> tuple[Optional[str], list[str]]:
        warnings: list[str] = []
        warnings.append("no suitable EBI coordinates result; attempting UniProt mapping")
        payload = {"from": "UniProtKB_AC-ID", "to": "Ensembl", "ids": uniprot_id}
        run = self.api.post(
            UNIPROT_IDMAP_RUN,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            data=payload,
        )
        if not isinstance(run.json_obj, dict):
            fallback, _ = self._fallback_uniprot_crossrefs(uniprot_id)
            if fallback:
                warnings.append("used UniProt cross-reference fallback mapping")
            return fallback, warnings
        job_id = run.json_obj.get("jobId") or run.json_obj.get("job_id")
        if not job_id:
            fallback, fallback_sources = self._fallback_uniprot_crossrefs(uniprot_id)
            if not fallback and fallback_sources:
                warnings.append(
                    "mapping API did not return jobId; available crossrefs include: "
                    f"{', '.join(sorted(set(fallback_sources)))}"
                )
            if fallback:
                warnings.append("used UniProt cross-reference fallback mapping")
                return fallback, warnings
            return None, warnings

        status = _poll_uniprot_status(self.api, str(job_id))
        if not status:
            warnings.append("UniProt mapping job did not complete")
            fallback, fallback_sources = self._fallback_uniprot_crossrefs(uniprot_id)
            if not fallback and fallback_sources:
                warnings.append(
                    "mapping job not completed; crossrefs available for "
                    f"{', '.join(sorted(set(fallback_sources)))}"
                )
            if fallback:
                warnings.append("used UniProt cross-reference fallback mapping")
                return fallback, warnings
            return None, warnings
        results = self.api.get(
            UNIPROT_IDMAP_RESULTS.format(job_id=job_id),
            headers={"Accept": "application/json"},
        )
        mapped_gene = self._extract_gene_from_mapping(results.json_obj)
        if mapped_gene:
            return mapped_gene, warnings
        fallback, fallback_sources = self._fallback_uniprot_crossrefs(uniprot_id)
        if not fallback and fallback_sources:
            warnings.append(
                "mapping results lacked Ensembl gene ID; crossrefs include "
                f"{', '.join(sorted(set(fallback_sources)))}"
            )
        if fallback:
            warnings.append("used UniProt cross-reference fallback mapping")
        return fallback, warnings

    def _fallback_uniprot_crossrefs(self, uniprot_id: str) -> tuple[Optional[str], list[str]]:
        try:
            response = self.api.get(
                UNIPROT_ENTRY_URL.format(accession=uniprot_id),
                headers={"Accept": "application/json"},
            )
        except ToolError:
            return None, []

        if not isinstance(response.json_obj, dict):
            return None, []

        refs = response.json_obj.get("uniProtKBCrossReferences", [])
        if not isinstance(refs, list):
            return None, []

        ensembl_sources: list[str] = []
        for ref in refs:
            if not isinstance(ref, dict):
                continue
            db_name = str(ref.get("database") or ref.get("dbDisplayName") or "").lower()
            if "ensembl" not in db_name:
                continue
            ensembl_sources.append(db_name)
            normalized = _extract_gene_from_ref_value(ref)
            if normalized:
                return normalized, ensembl_sources
        return None, ensembl_sources

    def _resolve_ncbi_gene(
        self,
        uniprot_id: str,
        uniprot_entry: Optional[dict[str, Any]] = None,
    ) -> tuple[Optional[dict[str, Any]], list[str]]:
        warnings: list[str] = []
        entry = uniprot_entry if isinstance(uniprot_entry, dict) else self._fetch_uniprot_entry(uniprot_id)
        if not isinstance(entry, dict):
            return None, ["unable to read UniProt entry for NCBI fallback"]

        organism = entry.get("organism") or {}
        taxid = _to_int(organism.get("taxonId"))
        organism_name = self._coerce_str(organism.get("scientificName"))
        gene_aliases = self._collect_ncbi_gene_aliases(entry)
        if not gene_aliases:
            return None, ["no gene alias found for NCBI fallback"]

        ncbi_accessions = self._collect_ncbi_nucleotide_accessions(entry)
        if ncbi_accessions:
            accession_preview = self._format_hint_list(ncbi_accessions, limit=12)
            warnings.append(
                f"NCBI fallback: collected {len(ncbi_accessions)} RefSeq/EMBL accession hint(s): {accession_preview}"
            )

        summaries = self._collect_ncbi_gene_summaries(gene_aliases, taxid, organism_name)
        if not summaries:
            return None, ["NCBI fallback: no NCBI gene record found for candidate identifiers"]

        chosen = self._choose_best_ncbi_gene_summary(summaries, gene_aliases, ncbi_accessions)
        if not chosen:
            return None, ["NCBI fallback: gene records lacked usable genomic location"]

        ncbi_coordinates = self._ncbi_summary_to_coordinates(chosen, gene_aliases, organism_name)
        if not ncbi_coordinates:
            return None, ["NCBI fallback: failed to extract genomic coordinates from chosen NCBI summary"]
        return ncbi_coordinates, warnings + ["NCBI fallback: resolved via NCBI gene summary"]

    def _fetch_uniprot_entry(self, uniprot_id: str) -> Optional[dict[str, Any]]:
        response = self.api.get(
            UNIPROT_ENTRY_URL.format(accession=uniprot_id),
            headers={"Accept": "application/json"},
        )
        if not isinstance(response.json_obj, dict):
            return None
        return response.json_obj

    def _collect_ncbi_gene_aliases(self, entry: dict[str, Any]) -> list[str]:
        aliases: list[str] = []

        for gene in entry.get("genes", []) or []:
            if not isinstance(gene, dict):
                continue

            gene_name = self._coerce_str(gene.get("geneName", {}).get("value"))
            if gene_name and gene_name not in aliases:
                aliases.append(gene_name)

            for key in ("orderedLocusNames", "orfNames"):
                raw = gene.get(key)
                if not isinstance(raw, list):
                    raw = [raw]
                for item in raw:
                    if isinstance(item, dict):
                        value = self._coerce_str(item.get("value"))
                    else:
                        value = self._coerce_str(item)
                    if value and value not in aliases:
                        aliases.append(value)

            for key in ("synonyms", "name"):
                raw = gene.get(key)
                if raw is None:
                    continue
                if isinstance(raw, list):
                    iterable = raw
                else:
                    iterable = [raw]
                for item in iterable:
                    if isinstance(item, dict):
                        value = self._coerce_str(item.get("value"))
                    else:
                        value = self._coerce_str(item)
                    if value and value not in aliases:
                        aliases.append(value)

        for ref in entry.get("uniProtKBCrossReferences", []) or []:
            if not isinstance(ref, dict):
                continue
            db_name = str(ref.get("database") or ref.get("dbDisplayName") or "").lower()
            if "ensemblbacteria" not in db_name:
                continue
            properties = self._coerce_properties(ref)
            for key in ("geneid", "genesymbol", "proteinid"):
                value = self._coerce_str(properties.get(key))
                if value and value not in aliases:
                    aliases.append(value)

        return aliases

    def _collect_ncbi_nucleotide_accessions(self, entry: dict[str, Any]) -> list[str]:
        accessions: list[str] = []
        for ref in entry.get("uniProtKBCrossReferences", []) or []:
            if not isinstance(ref, dict):
                continue
            db_name = str(ref.get("database") or ref.get("dbDisplayName") or "").lower()
            props = self._coerce_properties(ref)
            if "refseq" in db_name and props.get("nucleotidesequenceid"):
                accessions.append(props["nucleotidesequenceid"])
            elif db_name == "embl":
                molecule_type = (props.get("moleculetype") or "").lower()
                accession = self._coerce_str(ref.get("id"))
                if accession and accession not in accessions and "genomic_dna" in molecule_type.replace(" ", "_"):
                    accessions.append(accession)
        return accessions

    def _collect_ncbi_gene_summaries(
        self,
        aliases: list[str],
        taxid: Optional[int],
        organism_name: str,
    ) -> list[dict[str, Any]]:
        gene_ids: list[str] = []
        for alias in aliases:
            for term in self._build_ncbi_gene_terms(alias, taxid, organism_name):
                try:
                    ids = self._ncbi_esearch_gene_ids(term)
                except ToolError:
                    continue
                for gene_id in ids:
                    if gene_id not in gene_ids:
                        gene_ids.append(gene_id)
                if ids:
                    break
            if gene_ids:
                break

        if not gene_ids:
            return []

        try:
            return self._ncbi_gene_summaries(gene_ids)
        except ToolError:
            return []

    def _build_ncbi_gene_terms(self, alias: str, taxid: Optional[int], organism_name: str) -> list[str]:
        token = self._coerce_str(alias)
        if not token:
            return []
        terms: list[str] = []
        if taxid is not None:
            terms.append(f"{token}[Gene Name] AND {taxid}[Taxonomy ID]")
            terms.append(f"{token}[All Fields] AND {taxid}[Taxonomy ID]")
        if organism_name:
            genus_species = self._extract_genus_species(organism_name)
            if genus_species:
                terms.append(f'{token}[Gene Name] AND "{genus_species}"[Organism]')
        terms.append(f"{token}[Gene Name]")
        return terms

    def _ncbi_esearch_gene_ids(self, term: str) -> list[str]:
        response = self.api.get(
            NCBI_ESEARCH,
            headers={"Accept": "application/json"},
            params={"db": "gene", "term": term, "retmode": "json", "retmax": 50},
        )
        if not isinstance(response.json_obj, dict):
            return []
        result = response.json_obj.get("esearchresult")
        if not isinstance(result, dict):
            return []
        raw_ids = result.get("idlist", [])
        if not isinstance(raw_ids, list):
            return []
        ids: list[str] = []
        for raw_id in raw_ids:
            id_text = self._coerce_str(raw_id)
            if id_text and id_text not in ids:
                ids.append(id_text)
        return ids

    def _ncbi_gene_summaries(self, gene_ids: list[str]) -> list[dict[str, Any]]:
        response = self.api.get(
            NCBI_ESUMMARY,
            headers={"Accept": "application/json"},
            params={"db": "gene", "id": ",".join(gene_ids), "retmode": "json"},
        )
        if not isinstance(response.json_obj, dict):
            return []
        result = response.json_obj.get("result")
        if not isinstance(result, dict):
            return []
        uid_list = result.get("uids")
        if not isinstance(uid_list, list):
            return []
        summaries: list[dict[str, Any]] = []
        for uid in uid_list:
            record = result.get(uid)
            if isinstance(record, dict):
                record.setdefault("uid", uid)
                summaries.append(record)
        return summaries

    def _choose_best_ncbi_gene_summary(
        self,
        summaries: list[dict[str, Any]],
        aliases: list[str],
        accessions: list[str],
    ) -> Optional[dict[str, Any]]:
        target_accessions = {self._coerce_str(acc).upper() for acc in accessions if self._coerce_str(acc)}
        alias_set = {self._coerce_str(alias).lower() for alias in aliases if self._coerce_str(alias)}

        for summary in summaries:
            genomic_info = self._coerce_ncbi_gene_genomic_info(summary)
            if not genomic_info:
                continue
            if self._coerce_str(genomic_info.get("chraccver")).upper() in target_accessions:
                return summary

        for summary in summaries:
            genomic_info = self._coerce_ncbi_gene_genomic_info(summary)
            if not genomic_info:
                continue
            summary_aliases = self._coerce_ncbi_summary_aliases(summary)
            if summary_aliases.intersection(alias_set):
                return summary

        for summary in summaries:
            if self._coerce_ncbi_gene_genomic_info(summary):
                return summary
        return None

    def _coerce_ncbi_summary_aliases(self, summary: dict[str, Any]) -> set[str]:
        aliases: set[str] = set()
        for key in ("name", "nomenclaturesymbol", "nomenclaturename", "otherdesignations"):
            value = self._coerce_str(summary.get(key))
            if value:
                aliases.add(value.lower())
        other_aliases = self._coerce_str(summary.get("otheraliases"))
        if other_aliases:
            for item in other_aliases.split(","):
                alias = self._coerce_str(item)
                if alias:
                    aliases.add(alias.lower())
        return aliases

    def _coerce_ncbi_gene_genomic_info(self, summary: dict[str, Any]) -> Optional[dict[str, Any]]:
        genomicinfo = summary.get("genomicinfo")
        if not isinstance(genomicinfo, list):
            return None
        for item in genomicinfo:
            if not isinstance(item, dict):
                continue
            chr_start = _to_int(item.get("chrstart"))
            chr_stop = _to_int(item.get("chrstop"))
            chr_accver = self._coerce_str(item.get("chraccver"))
            if chr_start is not None and chr_stop is not None and chr_accver:
                return item
        return None

    def _ncbi_summary_to_coordinates(
        self,
        summary: dict[str, Any],
        aliases: list[str],
        organism_name: str,
    ) -> Optional[dict[str, Any]]:
        genomic_info = self._coerce_ncbi_gene_genomic_info(summary)
        if not genomic_info:
            return None
        raw_start = _to_int(genomic_info.get("chrstart"))
        raw_end = _to_int(genomic_info.get("chrstop"))
        if raw_start is None or raw_end is None:
            return None

        start = min(raw_start, raw_end)
        end = max(raw_start, raw_end)
        if raw_start <= raw_end:
            strand = 1
        else:
            strand = -1
        organism = summary.get("organism") or {}
        ncbi_accession = self._coerce_str(genomic_info.get("chraccver"))
        species_name = self._coerce_str(organism.get("scientificname")) or organism_name or "unknown"

        return {
            "ensembl_gene_id": self._coerce_str(summary.get("name"))
            or (self._coerce_str(aliases[0]) if aliases else "ncbi_gene"),
            "ncbi_accession": ncbi_accession or "unknown",
            "species": species_name,
            "assembly_name": f"NCBI {species_name}",
            "seq_region_name": ncbi_accession or "unknown",
            "gene_start_1based": start,
            "gene_end_1based": end,
            "strand": strand,
            "display_name": self._coerce_str(summary.get("name")) or self._coerce_str(summary.get("nomenclaturesymbol")),
            "taxid": _to_int(organism.get("taxid")),
        }

    def _coerce_str(self, value: Any) -> str:
        if value is None:
            return ""
        if isinstance(value, str):
            return value.strip()
        return str(value).strip()

    def _format_hint_list(self, values: list[str], limit: int = 12) -> str:
        unique: list[str] = []
        for raw in values:
            value = self._coerce_str(raw).upper()
            if value and value not in unique:
                unique.append(value)
        if len(unique) <= limit:
            return ", ".join(unique)
        extra = len(unique) - limit
        head = unique[:limit]
        return ", ".join(head) + f", ... (+{extra} more)"

    def _coerce_properties(self, entry_ref: dict[str, Any]) -> dict[str, str]:
        properties: dict[str, str] = {}
        for prop in entry_ref.get("properties", []) or []:
            if not isinstance(prop, dict):
                continue
            key = self._coerce_str(prop.get("key")).lower().replace(" ", "")
            value = self._coerce_str(prop.get("value"))
            if key and value:
                properties[key] = value
        return properties

    def _extract_genus_species(self, organism_name: str) -> str:
        base = organism_name.split("(")[0].strip()
        parts = base.split()
        if len(parts) >= 2:
            return f"{parts[0]} {parts[1]}"
        return ""

    def _lookup_ensembl_gene(self, ensembl_gene_id: str, _depth: int = 0) -> Optional[dict[str, Any]]:
        ensembl_gene_id = _normalize_ensembl_gene_id(ensembl_gene_id) or ensembl_gene_id
        resp = self.api.get(
            ENSEMBL_LOOKUP.format(ensembl_id=ensembl_gene_id),
            headers={"Accept": "application/json"},
            params={"expand": 0},
        )
        if not isinstance(resp.json_obj, dict):
            return None
        data = resp.json_obj
        if _depth > 3:
            return None
        object_type = str(data.get("object_type") or "").lower()
        if object_type and object_type != "gene":
            parent = data.get("Parent") or data.get("parent")
            if isinstance(parent, list):
                for item in parent:
                    if isinstance(item, str):
                        parent_id = _normalize_ensembl_gene_id(item)
                        if parent_id and parent_id != ensembl_gene_id:
                            return self._lookup_ensembl_gene(parent_id, _depth + 1)
                    if isinstance(item, dict):
                        parent_id = _normalize_ensembl_gene_id(item.get("id"))
                        if parent_id and parent_id != ensembl_gene_id:
                            return self._lookup_ensembl_gene(parent_id, _depth + 1)
            elif isinstance(parent, str):
                parent_id = _normalize_ensembl_gene_id(parent)
                if parent_id and parent_id != ensembl_gene_id:
                    return self._lookup_ensembl_gene(parent_id, _depth + 1)
            elif isinstance(parent, dict):
                parent_id = _normalize_ensembl_gene_id(parent.get("id"))
                if parent_id and parent_id != ensembl_gene_id:
                    return self._lookup_ensembl_gene(parent_id, _depth + 1)
        species = data.get("species")
        if isinstance(species, dict):
            species = species.get("name") or species.get("display_name") or species.get("scientific_name")
        return {
            "ensembl_gene_id": ensembl_gene_id,
            "species": str(species or "unknown"),
            "assembly_name": str(data.get("assembly_name", "unknown")),
            "seq_region_name": str(data.get("seq_region_name", data.get("seq_region", ""))),
            "gene_start_1based": _to_int(data.get("start")),
            "gene_end_1based": _to_int(data.get("end")),
            "strand": _to_int(data.get("strand"), default=1),
            "display_name": data.get("display_name") or data.get("external_name"),
            "taxid": _to_int(data.get("taxonomy_id") or data.get("taxid")),
        }

    def _extract_gene_from_mapping(self, mapping_payload: Any) -> Optional[str]:
        if not isinstance(mapping_payload, dict):
            return None
        entries = mapping_payload.get("results") or []
        if isinstance(mapping_payload.get("to"), str) and _normalize_ensembl_gene_id(mapping_payload.get("to")):
            return _normalize_ensembl_gene_id(mapping_payload.get("to"))
        if not isinstance(entries, list):
            entries = []
        for item in entries:
            to = self._normalize_mapping_value(item.get("to"))
            if to:
                return to
            to = self._normalize_mapping_value(item.get("toPrimaryAccession"))
            if to:
                return to
            to = self._normalize_mapping_value(item.get("to_id"))
            if to:
                return to
            to = self._normalize_mapping_value(item.get("toSecondary"))
            if to:
                return to
        for item in entries:
            for key in (item.get("to"), item.get("toSecondary"), item.get("to_id")):
                to = self._normalize_mapping_value(key)
                if to:
                    return to
        return None

    def _normalize_mapping_value(self, raw: Any) -> Optional[str]:
        return _normalize_ensembl_gene_id(raw)


def _to_int(value: Any, default: Optional[int] = None) -> Optional[int]:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _poll_uniprot_status(api: ApiClient, job_id: str, max_attempts: int = 20, delay_seconds: float = 1.0) -> bool:
    import time

    for _ in range(max_attempts):
        status = api.get(UNIPROT_IDMAP_STATUS.format(job_id=job_id), headers={"Accept": "application/json"})
        data = status.json_obj or {}
        if isinstance(data, dict):
            if "results" in data:
                results = data.get("results")
                if results:
                    return True
                failed_ids = data.get("failedIds")
                if failed_ids:
                    return False
            if data.get("jobStatus") in {"FINISHED", "COMPLETED", "FINISHED_WITH_WARNINGS"}:
                return bool(data.get("results"))
            if data.get("jobStatus") == "ERROR":
                return False
        time.sleep(delay_seconds)
    return False
