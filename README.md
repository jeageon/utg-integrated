# UTG (UniPcrTemplate)

`UTG`는 UniProtKB accession 하나로 유전자 주변 gDNA를 조회해 PCR 설계 회피 구간(negative feature)을 계산하고,
GenBank 파일(`.negfeatures.gb`)을 생성하는 CLI 프로그램입니다.

## 설치

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## 실행

```bash
python -m src.main P12345 --outdir data/output --flank 10000
```

주요 옵션:
- `--flank`: 유전자 좌우 확장 bp(기본 10000)
- `--flank-mode`: `genomic` 또는 `strand_relative`
- `--features`: `repeat,simple,variation,structural_variation,extreme_gc,homopolymer,ambiguous`
- `--mask`: `none|soft|hard` (Ensembl repeat mask)
- `--maf-threshold`: 변이 MAF 임계값
- `--gc-window`, `--gc-step`, `--gc-min`, `--gc-max`
- `--homopolymer-at`, `--homopolymer-gc`
- `--offline`: 캐시만 사용

출력 파일명:
`{UniProt}.{assembly}.{chr}_{extStart}_{extEnd}.negfeatures.gb`

동일 basename의 `...metadata.json`도 생성됩니다.

