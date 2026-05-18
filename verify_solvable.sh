#!/bin/bash

export GAP_DIR=$HOME/gap_experiment/gap-4.13.1
export PATH=$GAP_DIR:$PATH

echo "============================================================"
echo "Verifying 100 groups..."
echo "============================================================"

gap -q -b -T verify_solvable.g < /dev/null 2>&1

echo "============================================================"
echo "Done."
echo "============================================================"