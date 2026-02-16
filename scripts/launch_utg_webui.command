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

streamlit run src/webui.py
