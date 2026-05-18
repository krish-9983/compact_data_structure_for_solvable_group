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
orderG := Size(G);
Print("Benchmarking group: ", GROUP_EXPR, "\n");
Print("Order: ", Size(G), "\n");

iso  := IsomorphismPermGroup(G);
P    := Image(iso);
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
Print("Building index maps...\n");
gToP := ListWithIdenticalEntries(Length(els), 0);
pToG := ListWithIdenticalEntries(Length(els), 0);
for i in [1 .. Length(els)] do
    pIdx    := Position(elms, Image(iso, els[i])) - 1;  # 0-based perm index
    gToP[i] := pIdx;
    pToG[pIdx + 1] := i;
od;
Print("Index maps built.\n");

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
    while true do
        line := ReadLine(file);
        if line = fail then break; fi;
        line := Filtered(line, c -> c <> '\r' and c <> '\n');
        if Length(line) < 1 then continue; fi;
        tag  := line[1];
        # Parse all integers from the rest of the line.
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
        # Build the instruction record based on opcode character.
        if   tag = 'C' then Add(ops, rec(type := "const", id := nums[1]));
        elif tag = 'I' then Add(ops, rec(type := "inv",   a  := nums[1]));
        elif tag = 'M' then Add(ops, rec(type := "mul",   a  := nums[1], b := nums[2]));
        fi;
    od;
    CloseStream(file);
    return ops;
end;

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
EvaluateSLP := function(ops, leafValues)
    local values, i, op, inValA, inValB;
    values := ListWithIdenticalEntries(Length(ops), fail);
    for i in [1 .. Length(ops)] do
        op := ops[i];
        if   op.type = "const" then values[i] := leafValues[op.id + 1];
        elif op.type = "inv"   then values[i] := values[op.a + 1] ^ -1;
        else
            inValA    := values[op.a + 1];
            inValB    := values[op.b + 1];
            values[i] := inValA * inValB;    # GAP pc-group multiplication
        fi;
    od;
    return values[Length(values)];
end;

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
for i in [1 .. warmupReps] do
    EvaluateSLP(ops, leafValues);
od;

# Force a garbage collection before timing so GC pauses don't inflate the
# first timed rep. GASMAN("collect") is synchronous in GAP.
GASMAN("collect");

start := Runtime();   # GAP's millisecond wall clock
for i in [1 .. repetitions] do
    result := EvaluateSLP(ops, leafValues);
od;
finish := Runtime();

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
while Length(fracStr) > 0 and fracStr[Length(fracStr)] = '0' do
    fracStr := fracStr{[1..Length(fracStr)-1]};
od;

if Length(fracStr) = 0 then timeStr := String(intPart);
else timeStr := Concatenation(String(intPart), ".", fracStr);
fi;

Print("\n--- GAP Results ---\n");
Print("Result index in P (0-based): ", resultIndexInP, "\n");
Print("Total time (ms): ", totalTime, "\n");
Print("Time per repetition (ms): ", timePerRep, "\n");

# =============================================================================
# Append CSV row

# =============================================================================
filepath  := "benchmark_results.csv";

safeGroup := ReplacedString(GROUP_EXPR, ",", "_");
safeGroup := ReplacedString(safeGroup,  "(", "");
safeGroup := ReplacedString(safeGroup,  ")", "");

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

AppendTo(filepath, line);
Print("CSV appended successfully.\n");
