from __future__ import annotations

from ..config import (
    ENSEMBL_SEQUENCE_ID,
    ENSEMBL_SEQUENCE_REGION,
    ENSEMBL_SEQUENCE_MAX_BP,
    ENSEMBL_SEQUENCE_SAFETY_BP,
)
from ..models.data_schemas import GenomicCoordinates
from ..utils.api_client import ApiClient
from ..utils.coord_utils import clamp_region_length, build_region_string
from ..utils.exceptions import SequenceLengthMismatchError


class SequenceFetcher:
    def __init__(self, api_client: ApiClient) -> None:
        self.api = api_client

    def fetch(
        self,
        coordinates: GenomicCoordinates,
        mask: str = "soft",
        strict_region: bool = False,
    ) -> tuple[str, list[str]]:
        del strict_region
        warnings: list[str] = []
        region_len = coordinates.ext_end_1based - coordinates.ext_start_1based + 1
        region_start = coordinates.ext_start_1based
        region_end = coordinates.ext_end_1based
        if region_len > ENSEMBL_SEQUENCE_MAX_BP:
            region_start, region_end, clamped = clamp_region_length(
                region_start,
                region_end,
                ENSEMBL_SEQUENCE_SAFETY_BP,
            )
            if clamped:
                warnings.append(
                    "sequence span exceeded API limit; flank automatically reduced for safety"
                )
                region_len = region_end - region_start + 1

        seq, used = self._fetch_region(
            species=coordinates.species,
            chr_name=coordinates.seq_region_name,
            start_1based=region_start,
            end_1based=region_end,
            strand=coordinates.strand,
            mask=mask,
        )
        expected = region_end - region_start + 1
        if len(seq) != expected:
            raise SequenceLengthMismatchError(
                f"Expected {expected} nt but received {len(seq)} nt from Ensembl"
            )
        if coordinates.ext_start_1based != region_start or coordinates.ext_end_1based != region_end:
            used.append("sequence region was internally clamped by max request size")
        return seq.upper(), used

    def _fetch_region(
        self,
        species: str,
        chr_name: str,
        start_1based: int,
        end_1based: int,
        strand: int,
        mask: str,
    ) -> tuple[str, list[str]]:
        region = build_region_string(chr_name, start_1based, end_1based, strand)
        params = {}
        if mask and mask != "none":
            params["mask"] = mask
        resp = self.api.get(
            ENSEMBL_SEQUENCE_REGION.format(species=species, region=region),
            headers={"Content-Type": "text/plain"},
            params=params,
        )
        return self._normalize_fasta_or_plain(resp.text or ""), []

    def _normalize_fasta_or_plain(self, text: str) -> tuple[str, list[str]]:
        lines = [line.strip() for line in text.splitlines() if line.strip()]
        if not lines:
            return "", []
        if lines[0].startswith(">"):
            body = "".join(lines[1:])
            return body, ["input treated as FASTA; header removed"]
        return "".join(lines), []

