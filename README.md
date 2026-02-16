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

### Web UI 실행

```bash
streamlit run src/webui.py
```

브라우저에서 UniProt ID와 파이프라인 옵션을 입력하면 실행 결과(.gb, metadata.json)를 바로 다운로드할 수 있습니다.

### 바탕화면 바로가기(Quick Launch)

- 실행 스크립트:
  - `scripts/launch_utg_webui.command` (Web UI)
  - `scripts/launch_utg_cli.command` (CLI)
- 바탕화면에 바로가기 생성 예시:
  ```bash
  cp scripts/launch_utg_webui.command ~/Desktop/UTG-WebUI.command
  cp scripts/launch_utg_cli.command ~/Desktop/UTG-CLI.command
  ```
- `.command` 파일은 더블클릭 시 터미널에서 실행됩니다.

※ 현재 환경에서는 위 파일을 macOS 바탕화면에 복사해 두어 바로 실행 가능합니다.

주의:
- 웹 UI는 실행한 현재 머신에서 API 호출을 수행합니다.
- 출력 디렉터리는 웹앱 컨테이너의 파일시스템 경로를 사용합니다.

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
