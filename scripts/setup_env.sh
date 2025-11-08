#!/usr/bin/env bash
set -euo pipefail

ENV_FILE="env/environment.yml"
ENV_NAME=$(grep '^name:' "$ENV_FILE" | awk '{print $2}')

# Detect micromamba or mamba
if command -v micromamba >/dev/null 2>&1; then
    TOOL="micromamba"
elif command -v mamba >/dev/null 2>&1; then
    TOOL="mamba"
elif command -v conda >/dev/null 2>&1; then
    TOOL="conda"
else
    echo "ERROR: micromamba, mamba, or conda not found. Please install micromamba first:"
    echo "curl micro.mamba.pm/install.sh | bash"
    exit 1
fi

echo "🔧 Creating or updating environment using $TOOL..."
$TOOL env create -f "$ENV_FILE" || $TOOL env update -f "$ENV_FILE"

# Activate environment temporarily to install R packages from GitHub
echo "📦 Installing additional R packages (TESS3r)..."
if [ "$TOOL" = "micromamba" ]; then
    eval "$(micromamba shell hook -s bash)"
fi
$TOOL activate "$ENV_NAME"

Rscript -e 'if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools", repos="https://cloud.r-project.org")'
Rscript -e 'devtools::install_github("bcm-uga/TESS3_encho_sen", dependencies=TRUE, upgrade="never")'

echo "✅ Environment setup complete!"
echo "To activate the environment, run:"
echo "  micromamba activate $ENV_NAME"
echo
echo "Once activated, test with:"
echo "  malaria-pipeline --help"
