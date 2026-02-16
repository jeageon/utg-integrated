from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field


class GenomicCoordinates(BaseModel):
    uniprot_id: str
    ensembl_gene_id: str
    species: str
    assembly_name: str
    seq_region_name: str
    gene_start_1based: int
    gene_end_1based: int
    strand: int
    display_name: str | None = None
    taxid: int | None = None

    ext_start_1based: int
    ext_end_1based: int


class NegativeFeature(BaseModel):
    feature_type: str
    start: int
    end: int
    description: str
    source: str
    score: float | None = None
    strand: int | None = None
    attributes: dict[str, Any] = Field(default_factory=dict)


class SequenceRecordBundle(BaseModel):
    coordinates: GenomicCoordinates
    full_sequence: str
    features: list[NegativeFeature]
    metadata: dict[str, Any]

