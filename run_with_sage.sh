#!/bin/bash
# Wrapper-Skript zum Ausführen von Python-Skripten mit SageMath

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_NAME="Hoffbauer Antiprimes.py"

# Führe das Skript mit SageMath-Python aus
sage -python "$SCRIPT_DIR/$SCRIPT_NAME" "$@"
