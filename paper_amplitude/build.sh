#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
pandoc manuscript.md -f markdown+tex_math_single_backslash -s --shift-heading-level-by=-1 --number-sections -o manuscript.tex
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex > build.log
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex >> build.log
