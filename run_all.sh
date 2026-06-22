#!/bin/bash
set -e

FIXED_REPS=50
WARMUP_REPS=5
GROW_L_MIN=1000
GROW_L_MAX=10000
N_SIZES=10

export PATH=$HOME/gap_experiment/gap-4.13.1:$PATH

SEEDS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)

GROW_L_SIZES=()
for i in $(seq 0 $((N_SIZES - 1))); do
    L_VAL=$(( GROW_L_MIN + i * (GROW_L_MAX - GROW_L_MIN) / (N_SIZES - 1) ))
    GROW_L_SIZES+=("$L_VAL")
done


echo "============================================================"
echo "FINAL EXPERIMENT - 100 groups, 20 fixed seeds per group"
echo "  FIXED_REPS = $FIXED_REPS  |  WARMUP_REPS = $WARMUP_REPS"
echo "  Seeds      = ${SEEDS[*]}"
echo "  L sizes    = ${GROW_L_SIZES[*]}"
echo "  Total runs = 100 x 20 x 10 = 20000"
echo "============================================================"

g++ -O3 -o preparedsa          preparedsa.cpp
g++ -O3 -o structure_generator structure_generator.cpp
g++ -O3 -o querycds             querycds.cpp
echo "Compiled OK"

rm -f benchmark_results_c++.csv benchmark_results.csv

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RESULTS_DIR="results_100groups_${TIMESTAMP}"
mkdir -p "$RESULTS_DIR"

echo "Group,Order,Mode,DL,BinSize,CayleySize" > group_info.csv

GROUP_NUM=0
SKIPPED=0
FAILED_GROUPS=""

while IFS= read -r GROUP; do
    [ -z "$GROUP" ] && continue
    GROUP_NUM=$(( GROUP_NUM + 1 ))
    echo ""
    echo "=============================="
    echo "GROUP $GROUP_NUM/100: $GROUP"
    echo "=============================="

    gap -q -b -c "GROUP_EXPR:=\"$GROUP\";" extractGroupInfo.g < /dev/null

    if [ ! -f group_order.txt ] || [ ! -s cayley_parseable.txt ]; then
        echo "  ERROR: GAP output missing - skipping"
        SKIPPED=$(( SKIPPED + 1 ))
        FAILED_GROUPS="$FAILED_GROUPS\n  $GROUP"
        continue
    fi

    GROUP_ORDER=$(cat group_order.txt)
    ORDER_MINUS_ONE=$(( GROUP_ORDER - 1 ))
    echo "  Order: $GROUP_ORDER"

    if ! ./preparedsa > /dev/null 2>&1 || [ ! -s precomputed.bin ]; then
        echo "  ERROR: preparedsa failed - skipping"
        SKIPPED=$(( SKIPPED + 1 ))
        FAILED_GROUPS="$FAILED_GROUPS\n  $GROUP"
        continue
    fi

    BIN_SIZE=$(wc -c < precomputed.bin)
    CAYLEY_SIZE=$(( GROUP_ORDER * GROUP_ORDER * 4 ))

    if grep -q "e_arr" precomputed.txt 2>/dev/null; then
        GROUP_MODE=5
    else
        GROUP_MODE=3
    fi

    GROUP_DL=$(gap -q -b -T \
        -c "G:=EvalString(\"$GROUP\"); Print(DerivedLength(G), \"\n\"); QUIT_GAP(0);" \
        < /dev/null 2>/dev/null | head -1 | tr -d '[:space:]')
    [ -z "$GROUP_DL" ] && GROUP_DL="NA"

    echo "\"$GROUP\",$GROUP_ORDER,$GROUP_MODE,$GROUP_DL,$BIN_SIZE,$CAYLEY_SIZE" >> group_info.csv
    echo "  Mode=$GROUP_MODE  DL=$GROUP_DL  BinSize=$BIN_SIZE"

    # Per-group seeds: group_num*1000 + base. Different SLPs per group,
    # fixed across all L values within the group (preserves prefix-lock).
    GROUP_SEEDS=()
    for BASE in "${SEEDS[@]}"; do
        GROUP_SEEDS+=( $(( GROUP_NUM * 1000 + BASE )) )
    done

    for L in "${GROW_L_SIZES[@]}"; do
        TREE_IDX=0
        for SEED in "${GROUP_SEEDS[@]}"; do
            TREE_IDX=$(( TREE_IDX + 1 ))

            ./structure_generator "$L" "$ORDER_MINUS_ONE" "$SEED" > /dev/null 2>&1

            if (( RANDOM % 2 == 0 )); then
                ./querycds "$GROUP" "$L" "$FIXED_REPS" "$SEED" "$WARMUP_REPS" > /dev/null 2>&1
                gap -q -b \
                    -c "GROUP_EXPR:=\"$GROUP\"; FIXED_REPS:=$FIXED_REPS; SLP_SEED:=$SEED; WARMUP_REPS:=$WARMUP_REPS;" \
                    tree_benchmark.g < /dev/null > /dev/null 2>&1
            else
                gap -q -b \
                    -c "GROUP_EXPR:=\"$GROUP\"; FIXED_REPS:=$FIXED_REPS; SLP_SEED:=$SEED; WARMUP_REPS:=$WARMUP_REPS;" \
                    tree_benchmark.g < /dev/null > /dev/null 2>&1
                ./querycds "$GROUP" "$L" "$FIXED_REPS" "$SEED" "$WARMUP_REPS" > /dev/null 2>&1
            fi

            echo "  L=$L  seed=$SEED  tree=$TREE_IDX/20 done"
        done
    done

    echo "  GROUP $GROUP_NUM complete."
done < groups.txt

cp benchmark_results_c++.csv "$RESULTS_DIR/" 2>/dev/null || true
cp benchmark_results.csv     "$RESULTS_DIR/" 2>/dev/null || true
cp group_info.csv            "$RESULTS_DIR/" 2>/dev/null || true

echo ""
echo "============================================================"
echo "EXPERIMENT COMPLETE"
echo "  Groups processed : $(( GROUP_NUM - SKIPPED )) / 100"
echo "  Skipped          : $SKIPPED"
if [ -n "$FAILED_GROUPS" ]; then
    echo "  Failed groups:"
    echo -e "$FAILED_GROUPS"
fi
echo "  Results in       : $RESULTS_DIR/"
echo "  Files:"
echo "    benchmark_results_c++.csv"
echo "    benchmark_results.csv"
echo "    group_info.csv"
echo "============================================================"
