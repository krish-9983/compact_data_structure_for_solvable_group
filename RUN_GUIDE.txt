==============================================================================
PROJECT RUN GUIDE
Compact Data Structures for Finite Solvable Groups
==============================================================================


WHAT THE CODE DOES
------------------
Three stages run one after the other:

  Stage 1 (GAP)   extractGroupInfo.g    ->  cayley_parseable.txt
  Stage 2 (C++)   preparedsa.cpp        ->  precomputed.bin
  Stage 3 (C++)   querycds.cpp          ->  benchmark_results_c++.csv
  Stage 3 (GAP)   tree_benchmark.g      ->  benchmark_results.csv


==============================================================================
PREREQUISITES
==============================================================================

  # GAP (any version 4.13 or later)
  sudo apt install gap

  # C++ compiler
  sudo apt install g++

  # Python (only needed for plots, not the benchmark)
  pip3 install pandas matplotlib openpyxl --break-system-packages


==============================================================================
QUICK TEST  --  one group
==============================================================================

Run these commands one by one. SmallGroup(64,267) is used as the test group.

  Step 1: Compile C++ files

    g++ -O3 -o preparedsa          preparedsa.cpp
    g++ -O3 -o structure_generator structure_generator.cpp
    g++ -O3 -o querycds             querycds.cpp

  Step 2: GAP extracts composition series and Cayley tables

    gap -q -b -c 'GROUP_EXPR:="SmallGroup(64,267)";' extractGroupInfo.g < /dev/null

  Step 3: C++ builds the CDS binary

    ./preparedsa

  Step 4: Generate a random expression tree (size 1000, seed 42)

    GROUP_ORDER=$(cat group_order.txt)
    ./structure_generator 1000 $((GROUP_ORDER - 1)) 42

  Step 5: Run the CDS program

    ./querycds "SmallGroup(64,267)" 1000 50 42 5
                ^group name         ^L  ^reps ^seed ^warmup

  Step 6: Run the GAP program on the same tree

    gap -q -b \
      -c 'GROUP_EXPR:="SmallGroup(64,267)"; FIXED_REPS:=50; SLP_SEED:=42; WARMUP_REPS:=5;' \
      tree_benchmark.g < /dev/null

After Step 5 and 6, check that the Result column matches in both CSVs.
That means the C++ answer and GAP answer are the same -- correctness confirmed.


==============================================================================
FULL 100-GROUP EXPERIMENT
==============================================================================

BEFORE RUNNING -- fix the GAP path
-----------------------------------
Open run_all.sh and verify_solvable.sh in a text editor.
Find this line near the top of both files:

    export PATH="$HOME/gap_experiment/gap-4.13.1:$PATH"

Change it to point to where GAP is actually installed on your machine.
If GAP is a system install (sudo apt install gap), just delete that line.


Step 1: Verify all 100 groups

    bash verify_solvable.sh

  This checks every group in groups.txt is solvable and non-abelian.
  It should print:  TOTAL PASS: 100, FAIL: 0
  If any group fails, fix groups.txt before continuing.


Step 2: Run the experiment

    chmod +x run_all.sh
    ./run_all.sh

  This will take 6-8 hours. To keep it running after you close the terminal:

    tmux new -s bench          (open a new session)
    ./run_all.sh | tee experiment.log
    Press Ctrl+B then D        (detach -- keeps running after window closes)
    tmux attach -t bench       (reattach later to check progress)

  When it finishes, a folder called results_YYYYMMDD_HHMMSS/ will appear with:
    benchmark_results_c++.csv  --  CDS timings
    benchmark_results.csv      --  GAP timings
    group_info.csv             --  mode, derived length, file sizes


Step 3: Generate plots

    bash plot_all.sh results_YYYYMMDD_HHMMSS

  Replace the folder name with your actual timestamped folder name.

  This produces 5 plots inside results_YYYYMMDD_HHMMSS/plots/:
    plot1a_fixedL_vs_order.png      Figure 6.1 in thesis
    plot2a_growingL_vs_order.png    Figure 6.2
    plot3a_grand_vs_order.png       Figure 6.3
    impact2_grouped_bar_mode.png    Figure 6.4
    extra4_space_comparison.png     Figure 6.5


==============================================================================
KNOWN BUGS TO FIX BEFORE SHARING THE CODE
==============================================================================

  1. run_all.sh  (progress counter line)
     Bug:  prints "GROUP 1/2000" instead of "GROUP 1/100"
     Cause: ${#SEEDS[@]}00 expands to 2000 (20 seeds x "00")
     Fix:  change ${#SEEDS[@]}00 to just 100

  2. run_all.sh  (message printed at the very end)
     Bug:  tells you to run "python3 plot_advisor.py" but that file does not exist
     Fix:  the correct command is:  bash plot_all.sh <results_dir>

  3. plot_all.sh  (line near the top)
     Bug:  "source $HOME/cds_env/bin/activate" will crash on a fresh machine
           if that virtual environment does not exist
     Fix:  either delete that line, or create the env first:
             python3 -m venv ~/cds_env
             source ~/cds_env/bin/activate
             pip install pandas matplotlib openpyxl

  4. run_all.sh and verify_solvable.sh  (GAP path)
     Bug:  GAP path is hardcoded to gap-4.13.1 but thesis used GAP 4.15.1
     Fix:  update the PATH line to match your installed version


==============================================================================
FILES REFERENCE
==============================================================================

  extractGroupInfo.g       GAP script. Computes composition series and writes
                           Cayley tables using a global element numbering.

  preparedsa.cpp           C++ preprocessor. Reads cayley_parseable.txt,
                           decides Case 1 or Case 2, builds all lookup tables,
                           saves to precomputed.bin as tagged binary blocks.

  querycds.cpp             C++ query program. Loads precomputed.bin, detects
                           the construction mode from the tags present, then
                           evaluates expression trees and records timing.

  structure_generator.cpp  Generates random expression trees. Same seed always
                           gives the same tree. Larger tree with same seed
                           extends the smaller one (prefix-locked).

  tree_benchmark.g         GAP-side benchmark. Reads the same structure.txt
                           as querycds. Result must match C++ for correctness.

  run_all.sh               Full 100-group experiment driver. Loops over all
                           groups x seeds x L values = 20,000 timed runs.

  verify_solvable.sh       Pre-run check. Confirms all 100 groups are solvable
                           and non-abelian before the main experiment starts.

  plot_all.sh              Generates all 5 thesis plots from the CSV results.

  export_excel.py          Called by plot_all.sh. Converts raw CSVs to Excel
                           files needed by plot_6.py.

  plot_6.py                Called by plot_all.sh. Generates the 3 line plots.

  groups.txt               List of 100 groups. One group expression per line.

==============================================================================
