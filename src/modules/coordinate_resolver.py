from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from ..config import (
    ENSEMBL_LOOKUP,
    EBI_COORDINATES_URL,
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


class CoordinateResolver:
    def __init__(self, api_client: ApiClient) -> None:
        self.api = api_client

    def resolve(
        self,
        uniprot_id: str,
        flank_bp: int,
        flank_mode: str = "genomic",
        assembly_preference: str = "auto",
        taxid_filter: int | None = None,
    ) -> ResolverResult:
        del assembly_preference
        warnings: list[str] = []
        ensembl_gene_id, fallback_warnings = self._resolve_ensembl_gene(uniprot_id)
        warnings.extend(fallback_warnings)
        if not ensembl_gene_id:
            raise NoMappingError(f"No Ensembl gene identifier found for {uniprot_id}")

        lookup = self._lookup_ensembl_gene(ensembl_gene_id)
        if not lookup:
            raise NoMappingError(f"lookup failed for {ensembl_gene_id}")
        if taxid_filter is not None and lookup.get("taxid") != taxid_filter:
            warnings.append(f"taxid mismatch: got {lookup.get('taxid')}, requested {taxid_filter}")

        gene_start = lookup["gene_start_1based"]
        gene_end = lookup["gene_end_1based"]
        strand = lookup["strand"]
        ext_start, ext_end = apply_flank(gene_start, gene_end, flank_bp, flank_mode, strand)

        return ResolverResult(
            coordinates=GenomicCoordinates(
                uniprot_id=uniprot_id,
                ensembl_gene_id=lookup["ensembl_gene_id"],
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

    def _resolve_ensembl_gene(self, uniprot_id: str) -> tuple[str | None, list[str]]:
        warnings: list[str] = []
        try:
            response = self.api.get(EBI_COORDINATES_URL.format(accession=uniprot_id), headers={"Accept": "application/json"})
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
        except ToolError:
            pass

        return self._fallback_uniprot_mapping(uniprot_id)

    def _extract_ensembl_gene_id(self, entry: dict[str, Any]) -> str | None:
        direct = entry.get("ensemblGeneId") or entry.get("ensembl_gene_id")
        if isinstance(direct, str) and direct:
            return direct
        for ref in entry.get("crossReferences", []) or []:
            db_name = str(ref.get("dbDisplayName", "")).lower()
            if "ensembl" not in db_name:
                continue
            xid = str(ref.get("id") or ref.get("accession") or "").strip()
            if xid.startswith("ENSG"):
                return xid
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

    def _fallback_uniprot_mapping(self, uniprot_id: str) -> tuple[str | None, list[str]]:
        warnings: list[str] = []
        payload = {"from": "UniProtKB_AC-ID", "to": "Ensembl", "ids": uniprot_id}
        run = self.api.post(UNIPROT_IDMAP_RUN, headers={"Content-Type": "application/json"}, json_payload=payload)
        if not isinstance(run.json_obj, dict):
            return None, warnings
        job_id = run.json_obj.get("jobId") or run.json_obj.get("job_id")
        if not job_id:
            return None, warnings

        status = _poll_uniprot_status(self.api, str(job_id))
        if not status:
            warnings.append("UniProt mapping job did not complete")
            return None, warnings
        results = self.api.get(UNIPROT_IDMAP_RESULTS.format(job_id=job_id), headers={"Content-Type": "application/json"})
        return self._extract_gene_from_mapping(results.json_obj), warnings

    def _lookup_ensembl_gene(self, ensembl_gene_id: str) -> dict[str, Any] | None:
        resp = self.api.get(
            ENSEMBL_LOOKUP.format(ensembl_id=ensembl_gene_id),
            headers={"Content-Type": "application/json"},
            params={"expand": 0},
        )
        if not isinstance(resp.json_obj, dict):
            return None
        data = resp.json_obj
        if data.get("object_type") and data.get("object_type") != "Gene":
            pass
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

    def _extract_gene_from_mapping(self, mapping_payload: Any) -> str | None:
        if not isinstance(mapping_payload, dict):
            return None
        entries = mapping_payload.get("results") or []
        if not isinstance(entries, list):
            entries = []
        for item in entries:
            to = (item.get("to") or item.get("toPrimaryAccession") or "").upper()
            if to.startswith("ENSG"):
                return to
        # fallback: ensembl gene via nested fields
        for item in entries:
            key = item.get("to") or item.get("toSecondary") or item.get("to_id")
            if isinstance(key, str) and key.upper().startswith("ENSG"):
                return key
        return None


def _to_int(value: Any, default: int | None = None) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _poll_uniprot_status(api: ApiClient, job_id: str, max_attempts: int = 8, delay_seconds: float = 1.5) -> bool:
    import time

    for _ in range(max_attempts):
        status = api.get(UNIPROT_IDMAP_STATUS.format(job_id=job_id), headers={"Content-Type": "application/json"})
        data = status.json_obj or {}
        if isinstance(data, dict):
            if data.get("jobStatus") in {"FINISHED", "COMPLETED", "FINISHED_WITH_WARNINGS"}:
                return True
            if data.get("jobStatus") == "ERROR":
                return False
        time.sleep(delay_seconds)
    return False

