#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="/Users/jg/Documents/UTG-integrated"
cd "$PROJECT_DIR"

if [ -d ".venv" ]; then
  source .venv/bin/activate
elif [ -d "venv" ]; then
  source venv/bin/activate
fi

if [ -z "${VIRTUAL_ENV:-}" ]; then
  if [ -x "$(command -v python3)" ]; then
    python3 -m pip install -r requirements.txt
  else
    python -m pip install -r requirements.txt
  fi
fi

read -p "UniProt ID: " UNIPROT_ID
if [ -z "$UNIPROT_ID" ]; then
  echo "UniProt ID가 비어 있습니다. 종료합니다."
  exit 1
fi

python -m src.main "$UNIPROT_ID"
