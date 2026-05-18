#!/bin/bash
# =============================================================================
# run_all.sh  —  Full 100-group benchmark driver
# =============================================================================
#
# PURPOSE
#   Automates the complete experiment: compiles C++, loops over all groups in
#   groups.txt, runs the three-stage pipeline (GAP analysis -> preprocessing
#   -> timing) for each group at 10 tree sizes and 20 seeds, and collects all
#   CSV output into a timestamped results folder.
#
# HOW TO RUN
#   chmod +x run_all.sh
#   ./run_all.sh

# PREREQUISITES
#   - GAP 4.13.1
#   - g++ with C++17 support
# =============================================================================

set -uo pipefail

# =============================================================================
# Experiment parameters
# =============================================================================

FIXED_REPS=50       # timed repetitions per (group, L, seed) 
WARMUP_REPS=5       # untimed warm-up reps to prime CPU caches and GAP dispatch
GROW_L_MIN=1000     # smallest SLP tree size
GROW_L_MAX=10000    # largest SLP tree size
N_SIZES=10          # number of L values (evenly spaced from min to max)

# =============================================================================
# GAP path If using a local GAP build 
# =============================================================================
export PATH="$HOME/gap_experiment/gap-4.13.1:$PATH"

# =============================================================================
# Fixed seeds
# All 20 seeds are generated ONCE here
# =============================================================================
SEEDS=(101 202 303 404 505 606 707 808 909 1010 1111 1212 1313 1414 1515 1616 1717 1818 1919 2020)

# =============================================================================
# Compute L values
# Evenly spaced from GROW_L_MIN to GROW_L_MAX with N_SIZES points.
# =============================================================================
GROW_L_SIZES=()
for i in $(seq 0 $((N_SIZES - 1))); do
    L_VAL=$(( GROW_L_MIN + i * (GROW_L_MAX - GROW_L_MIN) / (N_SIZES - 1) ))
    GROW_L_SIZES+=("$L_VAL")
done

echo "============================================================"
echo "BENCHMARK EXPERIMENT"
echo "  FIXED_REPS  = $FIXED_REPS  |  WARMUP_REPS = $WARMUP_REPS"
echo "  Seeds       = ${SEEDS[*]}"
echo "  L sizes     = ${GROW_L_SIZES[*]}"
echo "  Total runs  = 100 groups x ${#SEEDS[@]} seeds x $N_SIZES L values = $((100 * ${#SEEDS[@]} * N_SIZES))"
echo "============================================================"

# =============================================================================
# Compile C++ files
# =============================================================================
g++ -O3 -o preparedsa          preparedsa.cpp
g++ -O3 -o structure_generator structure_generator.cpp
g++ -O3 -o querycds             querycds.cpp
echo "Compiled OK"

# Clear any CSVs from a previous partial run so we don't mix old and new data.
rm -f benchmark_results_c++.csv benchmark_results.csv

# Create a timestamped results directory so multiple runs don't overwrite each other.
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RESULTS_DIR="results_100groups_${TIMESTAMP}"
mkdir -p "$RESULTS_DIR"

# group_info.csv collects metadata about each group.
# BinSize = actual size of precomputed.bin (empirical CDS space).
# These are used by plot_advisor.py for the space comparison plots.
echo "Group,Order,Mode,DL,BinSize,CayleySize" > group_info.csv

GROUP_NUM=0
SKIPPED=0
FAILED_GROUPS=""

# =============================================================================
# Main group loop
# =============================================================================
while IFS= read -r GROUP; do
    [ -z "$GROUP" ] && continue
    GROUP_NUM=$(( GROUP_NUM + 1 ))

    echo ""
    echo "=============================="
    echo "GROUP $GROUP_NUM/${#SEEDS[@]}00: $GROUP"
    echo "=============================="

    # --------------------------------------------------------------------------
    # Stage 1: GAP algebraic analysis
    # extractGroupInfo.g reads GROUP_EXPR, computes the composition series,
    # and writes cayley_parseable.txt + group_order.txt.
    # The "< /dev/null" prevents GAP from trying to read stdin and hanging.
    # --------------------------------------------------------------------------
    gap -q -b -c "GROUP_EXPR:=\"$GROUP\";" extractGroupInfo.g < /dev/null

    # If GAP failed silently (bad group ID, not solvable, memory issue, etc.),
    if [ ! -f group_order.txt ] || [ ! -s cayley_parseable.txt ]; then
        echo "  ERROR: GAP output missing — skipping $GROUP"
        SKIPPED=$(( SKIPPED + 1 ))
        FAILED_GROUPS="$FAILED_GROUPS\n  $GROUP"
        continue
    fi

    GROUP_ORDER=$(cat group_order.txt)
    ORDER_MINUS_ONE=$(( GROUP_ORDER - 1 ))  # max valid element index for structure_generator
    echo "  Order: $GROUP_ORDER"

    # --------------------------------------------------------------------------
    # Stage 2: C++ preprocessing
    # --------------------------------------------------------------------------
    if ! ./preparedsa > /dev/null 2>&1 || [ ! -s precomputed.bin ]; then
        echo "  ERROR: preparedsa failed — skipping $GROUP"
        SKIPPED=$(( SKIPPED + 1 ))
        FAILED_GROUPS="$FAILED_GROUPS\n  $GROUP"
        continue
    fi

    # Collect space metadata for this group.
    BIN_SIZE=$(wc -c < precomputed.bin)
    CAYLEY_SIZE=$(( GROUP_ORDER * GROUP_ORDER * 4 ))  # hypothetical full Cayley table bytes

    # Detect construction mode from precomputed.txt (the human-readable companion).
    # Mode 5 (Case 2, Theorem 4 nested inside Theorem 3) is signalled by the
    # presence of the "e_arr" tag. Mode 3 (Case 1, Theorem 3 alone) is the default.
    if grep -q "e_arr" precomputed.txt 2>/dev/null; then
        GROUP_MODE=5
    else
        GROUP_MODE=3
    fi

    # Use GAP for derived length
    # The -T flag prevents GAP from printing timing info that would corrupt the output.
    GROUP_DL=$(gap -q -b -T \
        -c "G:=EvalString(\"$GROUP\"); Print(DerivedLength(G), \"\n\"); QUIT_GAP(0);" \
        < /dev/null 2>/dev/null | head -1 | tr -d '[:space:]')
    [ -z "$GROUP_DL" ] && GROUP_DL="NA"

    echo "\"$GROUP\",$GROUP_ORDER,$GROUP_MODE,$GROUP_DL,$BIN_SIZE,$CAYLEY_SIZE" >> group_info.csv
    echo "  Mode=$GROUP_MODE  DL=$GROUP_DL  BinSize=${BIN_SIZE} bytes"

    # --------------------------------------------------------------------------
    # Stage 3: benchmark loop over all (L, seed) pairs
    # --------------------------------------------------------------------------
    for L in "${GROW_L_SIZES[@]}"; do
        TREE_IDX=0
        for SEED in "${SEEDS[@]}"; do
            TREE_IDX=$(( TREE_IDX + 1 ))

            # Generate the SLP for this (L, seed) pair.
            # Writes structure.txt in the current directory — both querycds and
            # tree_benchmark.g will read this same file next.
            ./structure_generator "$L" "$ORDER_MINUS_ONE" "$SEED" > /dev/null 2>&1

            # Randomise engine order for each (group, L, seed) triple.
            # If querycds always ran first, it would always benefit from caches
            # pre-warmed by process startup — and vice versa if GAP ran first.
            # $RANDOM is a bash built-in (0..32767); odd = CDS first, even = GAP first.
            if (( RANDOM % 2 == 0 )); then
                # CDS first, then GAP
                ./querycds "$GROUP" "$L" "$FIXED_REPS" "$SEED" "$WARMUP_REPS" > /dev/null 2>&1
                gap -q -b \
                    -c "GROUP_EXPR:=\"$GROUP\"; FIXED_REPS:=$FIXED_REPS; SLP_SEED:=$SEED; WARMUP_REPS:=$WARMUP_REPS;" \
                    tree_benchmark.g < /dev/null > /dev/null 2>&1
            else
                # GAP first, then CDS
                gap -q -b \
                    -c "GROUP_EXPR:=\"$GROUP\"; FIXED_REPS:=$FIXED_REPS; SLP_SEED:=$SEED; WARMUP_REPS:=$WARMUP_REPS;" \
                    tree_benchmark.g < /dev/null > /dev/null 2>&1
                ./querycds "$GROUP" "$L" "$FIXED_REPS" "$SEED" "$WARMUP_REPS" > /dev/null 2>&1
            fi

            echo "  L=$L  seed=$SEED  tree=$TREE_IDX/${#SEEDS[@]} done"
        done
    done

    echo "  GROUP $GROUP_NUM complete."
done < groups.txt   # groups.txt is in the same directory as this script

# =============================================================================
# Collect results into timestamped folder
# =============================================================================
cp benchmark_results_c++.csv "$RESULTS_DIR/" 2>/dev/null || true
cp benchmark_results.csv     "$RESULTS_DIR/" 2>/dev/null || true
cp group_info.csv            "$RESULTS_DIR/" 2>/dev/null || true

echo ""
echo "============================================================"
echo "EXPERIMENT COMPLETE"
echo "  Groups processed : $(( GROUP_NUM - SKIPPED )) / $GROUP_NUM"
echo "  Skipped          : $SKIPPED"
if [ -n "$FAILED_GROUPS" ]; then
    echo "  Failed groups:"
    echo -e "$FAILED_GROUPS"
fi
echo "  Results in       : $RESULTS_DIR/"
echo ""
echo "  Next step — generate all plots:"
echo "    python3 plot_advisor.py \\"
echo "      --cds  $RESULTS_DIR/benchmark_results_c++.csv \\"
echo "      --gap  $RESULTS_DIR/benchmark_results.csv \\"
echo "      --info $RESULTS_DIR/group_info.csv \\"
echo "      --out  $RESULTS_DIR/plots/"
echo "============================================================"
