#!/usr/bin/env bash
set -euo pipefail

ENV_FILE="env/environment.yml"
ENV_NAME=$(grep '^name:' "$ENV_FILE" | awk '{print $2}')

# Detect micromamba / mamba / conda
if command -v micromamba >/dev/null 2>&1; then
  TOOL="micromamba"
elif command -v mamba >/dev/null 2>&1; then
  TOOL="mamba"
elif command -v conda >/dev/null 2>&1; then
  TOOL="conda"
else
  echo "ERROR: micromamba, mamba, or conda not found."
  echo "Install micromamba first:"
  echo "  curl -L micro.mamba.pm/install.sh | bash"
  exit 1
fi

echo "🔧 Creating or updating environment using $TOOL..."
$TOOL env create -f "$ENV_FILE" || $TOOL env update -f "$ENV_FILE"

# Activate env for R/Python installs and compiler path
if [ "$TOOL" = "micromamba" ]; then
  eval "$(micromamba shell hook -s bash)"
fi
$TOOL activate "$ENV_NAME"

# ---- Install TESS3r (R) ----
echo "📦 Installing R package: TESS3r (bcm-uga/TESS3_encho_sen)"
Rscript -e 'if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools", repos="https://cloud.r-project.org")'
Rscript -e 'devtools::install_github("bcm-uga/TESS3_encho_sen", dependencies=TRUE, upgrade="never")'

# ---- Compile hmmIBD from source and install into env/bin ----
echo "🔨 Installing hmmIBD (from GitHub source)"
WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT
git -C "$WORKDIR" clone --depth 1 https://github.com/glipsnort/hmmIBD.git
cc -o "$WORKDIR/hmmIBD/hmmIBD" -O3 -Wall "$WORKDIR/hmmIBD/hmmIBD.c" -lm
INSTALL_BIN="${CONDA_PREFIX}/bin"
mkdir -p "$INSTALL_BIN"
cp "$WORKDIR/hmmIBD/hmmIBD" "$INSTALL_BIN/"
echo "✅ hmmIBD installed at: $INSTALL_BIN/hmmIBD"

# ---- Install your package so the CLI is available immediately ----
echo "📦 Installing malaria-popgen-toolkit (editable)"
pip install -e .

echo
echo "✅ Environment setup complete!"
echo "To activate, run:"
echo "  micromamba activate $ENV_NAME"
echo
echo "Smoke tests (optional):"
echo "  malaria-pipeline --help"
echo "  hmmIBD -h"
echo '  R -q -e "library(TESS3r); sessionInfo()"'

