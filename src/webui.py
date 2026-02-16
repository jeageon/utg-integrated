from __future__ import annotations

import json
from pathlib import Path

import streamlit as st

from src.config import DEFAULT_CACHE_TTL_HOURS, DEFAULT_FEATURES, DEFAULT_FLANK, DEFAULT_RETRIES, DEFAULT_TIMEOUT
from src.main import run_pipeline


st.set_page_config(page_title="UTG Web UI", page_icon="🧬", layout="wide")

st.title("UTG")
st.caption("UniPcrTemplate (UTG) - UniProt ID 기반 gDNA negative feature 추출")


with st.form("utg_form"):
    st.subheader("실행 설정")
    col1, col2, col3 = st.columns(3)

    with col1:
        uniprot_id = st.text_input("UniProt ID", placeholder="예: P12345")
    with col2:
        outdir = st.text_input("출력 디렉터리", value=str(Path("data/output")))
    with col3:
        flank = st.number_input("Flank (bp)", min_value=0, step=100, value=DEFAULT_FLANK)

    col4, col5, col6 = st.columns(3)
    with col4:
        flank_mode = st.selectbox("Flank mode", ["genomic", "strand_relative"], index=0)
    with col5:
        mask = st.selectbox("Mask", ["none", "soft", "hard"], index=1)
    with col6:
        assembly = st.selectbox("Assembly", ["auto", "GRCh38", "GRCh37"], index=0)

    feature_input = st.multiselect(
        "negative feature",
        options=[
            "repeat",
            "simple",
            "variation",
            "structural_variation",
            "extreme_gc",
            "homopolymer",
            "ambiguous",
        ],
        default=DEFAULT_FEATURES,
    )

    col7, col8, col9 = st.columns(3)
    with col7:
        maf_threshold = st.number_input("MAF threshold", min_value=0.0, max_value=1.0, value=0.01, step=0.01, format="%.4f")
    with col8:
        gc_window = st.number_input("GC window", min_value=1, value=50, step=1)
    with col9:
        gc_step = st.number_input("GC step", min_value=1, value=10, step=1)

    col10, col11, col12 = st.columns(3)
    with col10:
        gc_min = st.number_input("GC min", min_value=0.0, max_value=100.0, value=30.0, step=0.5)
    with col11:
        gc_max = st.number_input("GC max", min_value=0.0, max_value=100.0, value=70.0, step=0.5)
    with col12:
        homopolymer_at = st.number_input("Homopolymer A/T", min_value=2, value=5, step=1)

    homopolymer_gc = st.number_input("Homopolymer G/C", min_value=2, value=4, step=1)

    col13, col14, col15 = st.columns(3)
    with col13:
        timeout = st.number_input("Timeout", min_value=1.0, value=DEFAULT_TIMEOUT, step=0.5)
    with col14:
        retries = st.number_input("Retries", min_value=1, value=DEFAULT_RETRIES, step=1)
    with col15:
        cache_ttl_hours = st.number_input("Cache TTL (hour)", min_value=1, value=DEFAULT_CACHE_TTL_HOURS, step=1)

    offline = st.toggle("Offline mode", value=False, help="캐시만 사용")
    cache_mode = st.radio("Cache", options=["on", "off"], horizontal=True, index=0)
    write_metadata_json = st.toggle("metadata.json 저장", value=True)

    submitted = st.form_submit_button("UTG 실행", use_container_width=True)


if not submitted:
    st.info("좌측 폼을 작성하고 **UTG 실행** 버튼을 눌러주세요.")
    st.stop()


if not uniprot_id.strip():
    st.error("UniProt ID를 입력해주세요.")
    st.stop()


if gc_min > gc_max:
    st.error("GC min 값은 GC max 값보다 작아야 합니다.")
    st.stop()

with st.spinner("UniProt ID를 처리 중입니다..."):
    try:
        gb_path, metadata_path, summary = run_pipeline(
            uniprot_id=uniprot_id.strip(),
            outdir=Path(outdir),
            flank=int(flank),
            flank_mode=flank_mode,
            assembly=assembly,
            mask=mask,
            features=feature_input,
            maf_threshold=float(maf_threshold),
            gc_window=int(gc_window),
            gc_step=int(gc_step),
            gc_min=float(gc_min),
            gc_max=float(gc_max),
            homopolymer_at=int(homopolymer_at),
            homopolymer_gc=int(homopolymer_gc),
            timeout=float(timeout),
            retries=int(retries),
            cache=cache_mode,
            cache_ttl_hours=int(cache_ttl_hours),
            offline=offline,
            write_metadata_json=write_metadata_json,
        )
    except Exception as exc:
        st.error(f"실행 중 오류가 발생했습니다: {exc}")
        st.exception(exc)
        st.stop()

st.success("생성 완료")

left, right = st.columns(2)
with left:
    st.subheader("요약")
    st.write(f"UniProt ID: `{summary['uniprot_id']}`")
    st.write(f"GenBank: `{gb_path}`")
    st.write(f"총 feature: `{summary['n_features']}`")
    if summary["feature_counts"]:
        st.write("Feature count")
        st.bar_chart(summary["feature_counts"])

with right:
    st.subheader("좌표 정보")
    st.json(summary["coordinates"])

if summary.get("warnings"):
    with st.expander("Warnings"):
        for warning in summary["warnings"]:
            st.write(f"- {warning}")

gb_bytes = Path(summary["gb_path"]).read_bytes()
st.download_button(
    label="GenBank 파일 다운로드",
    data=gb_bytes,
    file_name=Path(summary["gb_path"]).name,
    mime="application/octet-stream",
    use_container_width=True,
)

if summary.get("metadata_path"):
    metadata_file = Path(summary["metadata_path"])
    metadata_text = metadata_file.read_text(encoding="utf-8")
    with st.expander("metadata.json 미리보기"):
        try:
            st.json(json.loads(metadata_text))
        except Exception:
            st.text(metadata_text)
    st.download_button(
        label="metadata.json 다운로드",
        data=metadata_text,
        file_name=metadata_file.name,
        mime="application/json",
        use_container_width=True,
    )

st.caption("결과 파일은 서버의 출력 폴더에 저장됩니다.")

