#!/usr/bin/env bash
set -euo pipefail

# Wechsle in das Verzeichnis, in dem dieses Skript liegt
cd "$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

# Entfernt typische LaTeX-Hilfsdateien im Skriptverzeichnis
rm -f -- \
  *.aux *.log *.out *.toc *.synctex.gz *.fdb_latexmk *.fls \
  *.bbl *.blg *.bcf *.run.xml *.nav *.snm *.vrb *.lof *.lot *.lol

echo "LaTeX-Hilfsdateien im Verzeichnis '$(pwd)' entfernt."
