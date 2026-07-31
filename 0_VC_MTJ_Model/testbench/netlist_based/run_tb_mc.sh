#!/usr/bin/env bash

set -e

# Run from the directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Remove results from the previous run
rm -rf psf_mc
rm -f spectre_mc.log

echo "Running VC-MTJ Monte Carlo Spectre testbench..."

spectre tb_mc.scs \
    +log spectre_mc.log \
    -raw ./psf_mc

echo "Simulation completed."
echo "Results: $SCRIPT_DIR/psf_mc"
echo "Opening ViVA..."

viva &
