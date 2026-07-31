#!/usr/bin/env bash

set -e

# Run from the directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Remove results from the previous run
rm -rf psf_typical
rm -f spectre_typical.log

echo "Running VC-MTJ typical-corner Spectre testbench..."

spectre tb_typical.scs \
    +log spectre_typical.log \
    -raw ./psf_typical

echo "Simulation completed."
echo "Results: $SCRIPT_DIR/psf_typical"
echo "Opening ViVA..."

viva &
