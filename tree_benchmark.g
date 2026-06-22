<<<<<<< HEAD
# 1. Construct Group
if not IsBound(GROUP_EXPR) then
    Error("Use: gap -q -c 'GROUP_EXPR:=\"<group>\";' tree_benchmark.g");
fi;

G := EvalString(GROUP_EXPR);
=======
# =============================================================================
# tree_benchmark.g  —  GAP-side SLP evaluator and benchmark timer
# =============================================================================
#
# PURPOSE
#   This is the GAP engine in the two-engine benchmark. It reads the same
#   structure.txt that querycds reads, evaluates the SLP using GAP's native
#   polycyclic group multiplication, and appends one timing  to
#   benchmark_results.csv.
#
# HOW TO CALL
#   gap -q -b \
#     -c "GROUP_EXPR:=\"SmallGroup(64,267)\"; FIXED_REPS:=50; SLP_SEED:=101; WARMUP_REPS:=5;" \
#     tree_benchmark.g
#
#   All four variables should be set. Defaults (50 reps, seed 0, 5 warmup).
#
# We use PC-GROUP FOR TIMING BUT PERM-GROUP FOR INDICES because : 
#   GAP's polycyclic collector (the pc-group multiplication) is what we want
#   to time — it's the "native group arithmetic" we're comparing against.
#   But the element indices in structure.txt use the permutation-group ordering
#   (same as C++), so we need pToG to convert those indices to pc-group elements
#   before evaluating the SLP, and gToP to convert the final result back to a
#   perm-group index for comparison with C++'s output.
# =============================================================================

# =============================================================================
# 1. Load and convert the group
# =============================================================================

if not IsBound(GROUP_EXPR) then
    Error("Use: gap -q -c 'GROUP_EXPR:=\"SmallGroup(96,227)\";' tree_benchmark.g");
fi;

G      := EvalString(GROUP_EXPR);
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
orderG := Size(G);
Print("Benchmarking group: ", GROUP_EXPR, "\n");
Print("Order: ", Size(G), "\n");

iso  := IsomorphismPermGroup(G);
P    := Image(iso);
<<<<<<< HEAD
els  := Elements(G);    # pc-group ordering (for timing)
elms := Elements(P);    # perm-group ordering (for index matching)

# Build gToP and pToG maps once (excluded from timing)
=======
els  := Elements(G);    # pc-group element list (for timing arithmetic)
elms := Elements(P);    # perm-group element list (for index matching with C++)

# =============================================================================
# 2. Build index translation maps (excluded from timing)
#
# gToP[i] : 0-based perm-group index of G's i-th element (1-based in GAP)
# pToG[j] : 1-based G-index of the element at perm-group index j (0-based)
#
# These maps let us convert a constant leaf ID from structure.txt (which is
# a perm-group 0-based index) into the corresponding pc-group element for
# evaluation, and then convert the final result back to a perm-group index
# so it matches what querycds writes to its CSV..
# =============================================================================
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
Print("Building index maps...\n");
gToP := ListWithIdenticalEntries(Length(els), 0);
pToG := ListWithIdenticalEntries(Length(els), 0);
for i in [1 .. Length(els)] do
<<<<<<< HEAD
    pIdx    := Position(elms, Image(iso, els[i])) - 1;  # 0-based P-index
=======
    pIdx    := Position(elms, Image(iso, els[i])) - 1;  # 0-based perm index
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
    gToP[i] := pIdx;
    pToG[pIdx + 1] := i;
od;
Print("Index maps built.\n");

<<<<<<< HEAD
# 2. Load structure.txt
ReadStructure := function(filename)
    local file, line, ops, tag, nums, p, n1;
    file := InputTextFile(filename);
    if file = fail then Error("Cannot open: ", filename); fi;
    ops  := [];
    line := ReadLine(file);
=======
# =============================================================================
# 3. SLP file reader
#
# Reads structure.txt and returns a list of instruction records.
# Each record has: type ("const"/"inv"/"mul"), and operand fields.
# =============================================================================
ReadStructure := function(filename)
    local file, line, ops, tag, nums, p, n1;
    file := InputTextFile(filename);
    if file = fail then Error("Cannot open structure file: ", filename); fi;
    ops  := [];
    line := ReadLine(file);   # skip the "L <size>" header line
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
    while true do
        line := ReadLine(file);
        if line = fail then break; fi;
        line := Filtered(line, c -> c <> '\r' and c <> '\n');
        if Length(line) < 1 then continue; fi;
        tag  := line[1];
<<<<<<< HEAD
=======
        # Parse all integers from the rest of the line.
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
        nums := [];
        p    := 1;
        while p <= Length(line) do
            if line[p] in "0123456789" then
                n1 := "";
                while p <= Length(line) and line[p] in "0123456789" do
                    n1 := Concatenation(n1, [line[p]]);
                    p  := p + 1;
                od;
                Add(nums, Int(n1));
            else
                p := p + 1;
            fi;
        od;
<<<<<<< HEAD
        if   tag = 'C' then
            Add(ops, rec(type := "const", id := nums[1]));
        elif tag = 'I' then
            Add(ops, rec(type := "inv",   a  := nums[1]));
        elif tag = 'M' then
            Add(ops, rec(type := "mul",   a  := nums[1], b := nums[2]));
=======
        # Build the instruction record based on opcode character.
        if   tag = 'C' then Add(ops, rec(type := "const", id := nums[1]));
        elif tag = 'I' then Add(ops, rec(type := "inv",   a  := nums[1]));
        elif tag = 'M' then Add(ops, rec(type := "mul",   a  := nums[1], b := nums[2]));
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
        fi;
    od;
    CloseStream(file);
    return ops;
end;

<<<<<<< HEAD
# 3. File path
STRUCTURE_FILE := "structure.txt";
ops := ReadStructure(STRUCTURE_FILE);

# 4. Build leaf values
# const IDs in structure.txt are P-indices (C++ ordering)
# pToG converts them to G-indices so SLP runs in pc-group
maxLeafId  := Maximum(List(Filtered(ops, o -> o.type = "const"), o -> o.id));
leafValues := List([0 .. maxLeafId], i -> els[pToG[i + 1]]);

# 5. Evaluate SLP
=======
# structure.txt is always read from the current working directory.
STRUCTURE_FILE := "structure.txt";
ops := ReadStructure(STRUCTURE_FILE);

# =============================================================================
# 4. Build leaf value array
#
# Constant leaf IDs in structure.txt are 0-based perm-group indices.
# We pre-build leafValues so that leafValues[id+1] is the pc-group element
# corresponding to perm-group index id. The +1 shifts from 0-based to GAP's
# 1-based list indexing.
# =============================================================================
maxLeafId  := Maximum(List(Filtered(ops, o -> o.type = "const"), o -> o.id));
leafValues := List([0 .. maxLeafId], i -> els[pToG[i + 1]]);

# =============================================================================
# 5. SLP evaluator
#
# Evaluates the SLP sequentially. Each instruction's result is stored in
# values[i] (1-based in GAP). Operand indices a and b in the instruction
# records are 0-based, so we always access values[a+1] and values[b+1].
# =============================================================================
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
EvaluateSLP := function(ops, leafValues)
    local values, i, op, inValA, inValB;
    values := ListWithIdenticalEntries(Length(ops), fail);
    for i in [1 .. Length(ops)] do
        op := ops[i];
<<<<<<< HEAD
        if op.type = "const" then
            values[i] := leafValues[op.id + 1];
        elif op.type = "inv" then
            values[i] := values[op.a + 1] ^ -1;
        else
            inValA    := values[op.a + 1];
            inValB    := values[op.b + 1];
            values[i] := inValA * inValB;
=======
        if   op.type = "const" then values[i] := leafValues[op.id + 1];
        elif op.type = "inv"   then values[i] := values[op.a + 1] ^ -1;
        else
            inValA    := values[op.a + 1];
            inValB    := values[op.b + 1];
            values[i] := inValA * inValB;    # GAP pc-group multiplication
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
        fi;
    od;
    return values[Length(values)];
end;

<<<<<<< HEAD
# 6. Benchmark - ONLY pc-group arithmetic timed
repetitions := 50; # default fixed repetitions
if IsBound(FIXED_REPS) then
    repetitions := FIXED_REPS;
fi;

slpSeed := 0; # default seed (for CSV tracking)
if IsBound(SLP_SEED) then
    slpSeed := SLP_SEED;
fi;

warmupReps := 5; # untimed warmup — primes GAP method cache
if IsBound(WARMUP_REPS) then
    warmupReps := WARMUP_REPS;
fi;

# Warmup — NOT timed.
# GAP interprets EvaluateSLP; the first call is slower because method
# caches are cold. Warmup reps run before the clock so timed reps
# reflect steady-state performance, not first-call compilation overhead.
=======
# =============================================================================
# 6. Benchmark parameters
# =============================================================================

# Fixed repetition count

slpSeed := 0;  # only written to CSV for tracking;
if IsBound(SLP_SEED) then slpSeed := SLP_SEED; fi;

warmupReps := 5;
if IsBound(WARMUP_REPS) then warmupReps := WARMUP_REPS; fi;

# =============================================================================
# 7. Timed benchmark
# =============================================================================

# Warm-up — NOT included in timing.
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
for i in [1 .. warmupReps] do
    EvaluateSLP(ops, leafValues);
od;

<<<<<<< HEAD
GASMAN("collect");
start := Runtime();
=======
# Force a garbage collection before timing so GC pauses don't inflate the
# first timed rep. GASMAN("collect") is synchronous in GAP.
GASMAN("collect");

start := Runtime();   # GAP's millisecond wall clock
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
for i in [1 .. repetitions] do
    result := EvaluateSLP(ops, leafValues);
od;
finish := Runtime();

<<<<<<< HEAD
# 7. Compute result index AFTER timing
resultIndexInG := Position(els, result) - 1;       # G ordering
resultIndexInP := gToP[resultIndexInG + 1];         # P ordering (matches C++)

# Verification
Print("Verification - result matches in P: ",
      elms[resultIndexInP + 1] = Image(iso, result), "\n");

totalTime  := finish - start;
# compute time per rep
timePerRep := Float(totalTime) / repetitions;

# multiply by 10^7 and round
scaled := Int(timePerRep * 10000000 + 0.5);

intPart := QuoInt(scaled, 10000000);
fracPart := scaled mod 10000000;

# build fractional string padded to 7 digits
fracStr := String(fracPart);
while Length(fracStr) < 7 do
    fracStr := Concatenation("0", fracStr);
od;

# remove trailing zeros
=======
# =============================================================================
# Compute result index AFTER timing
# =============================================================================
resultIndexInG := Position(els, result) - 1;        # 0-based G-index
resultIndexInP := gToP[resultIndexInG + 1];          # 0-based P-index (matches C++)

# Sanity check: the P-image of `result` should equal elms[resultIndexInP+1].
Print("Verification - result matches in P: ",
      elms[resultIndexInP + 1] = Image(iso, result), "\n");

# =============================================================================

totalTime  := finish - start;
timePerRep := Float(totalTime) / repetitions;

scaled   := Int(timePerRep * 10000000 + 0.5);
intPart  := QuoInt(scaled, 10000000);
fracPart := scaled mod 10000000;

fracStr := String(fracPart);
while Length(fracStr) < 7 do fracStr := Concatenation("0", fracStr); od;
# Strip trailing zeros so "5.3000000" becomes "5.3".
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
while Length(fracStr) > 0 and fracStr[Length(fracStr)] = '0' do
    fracStr := fracStr{[1..Length(fracStr)-1]};
od;

<<<<<<< HEAD
if Length(fracStr) = 0 then
    timeStr := String(intPart);
else
    timeStr := Concatenation(String(intPart), ".", fracStr);
=======
if Length(fracStr) = 0 then timeStr := String(intPart);
else timeStr := Concatenation(String(intPart), ".", fracStr);
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
fi;

Print("\n--- GAP Results ---\n");
Print("Result index in P (0-based): ", resultIndexInP, "\n");
Print("Total time (ms): ", totalTime, "\n");
Print("Time per repetition (ms): ", timePerRep, "\n");

<<<<<<< HEAD
# 8. Write CSV


filepath := "benchmark_results.csv";

safeGroup := ReplacedString(GROUP_EXPR, ",", "_");
safeGroup := ReplacedString(safeGroup, "(", "");
safeGroup := ReplacedString(safeGroup, ")", "");
=======
# =============================================================================
# Append CSV row

# =============================================================================
filepath  := "benchmark_results.csv";

safeGroup := ReplacedString(GROUP_EXPR, ",", "_");
safeGroup := ReplacedString(safeGroup,  "(", "");
safeGroup := ReplacedString(safeGroup,  ")", "");
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033

line := Concatenation(
    "GAP,", safeGroup, ",",
    String(orderG), ",",
    String(Length(ops)), ",",
    String(slpSeed), ",",
    String(repetitions), ",",
    String(resultIndexInP), ",",
    String(Int(totalTime)), ",",
    timeStr, "\n"
);

out := OutputTextFile(filepath, true);
SetPrintFormattingStatus(out, false);
WriteAll(out, line);
CloseStream(out);

Print("CSV appended successfully.\n");

# To run: Read("C:/Users/krish/Downloads/advanced_algorithm/goupCompactDS/tree_benchmark.g");
=======
Print("CSV appended successfully.\n");
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
