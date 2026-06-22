# 1. Construct Group
if not IsBound(GROUP_EXPR) then
    Error("Use: gap -q -c 'GROUP_EXPR:=\"<group>\";' tree_benchmark.g");
fi;

G := EvalString(GROUP_EXPR);
orderG := Size(G);
Print("Benchmarking group: ", GROUP_EXPR, "\n");
Print("Order: ", Size(G), "\n");

iso  := IsomorphismPermGroup(G);
P    := Image(iso);
els  := Elements(G);    # pc-group ordering (for timing)
elms := Elements(P);    # perm-group ordering (for index matching)

# Build gToP and pToG maps once (excluded from timing)
Print("Building index maps...\n");
gToP := ListWithIdenticalEntries(Length(els), 0);
pToG := ListWithIdenticalEntries(Length(els), 0);
for i in [1 .. Length(els)] do
    pIdx    := Position(elms, Image(iso, els[i])) - 1;  # 0-based P-index
    gToP[i] := pIdx;
    pToG[pIdx + 1] := i;
od;
Print("Index maps built.\n");

# 2. Load structure.txt
ReadStructure := function(filename)
    local file, line, ops, tag, nums, p, n1;
    file := InputTextFile(filename);
    if file = fail then Error("Cannot open: ", filename); fi;
    ops  := [];
    line := ReadLine(file);
    while true do
        line := ReadLine(file);
        if line = fail then break; fi;
        line := Filtered(line, c -> c <> '\r' and c <> '\n');
        if Length(line) < 1 then continue; fi;
        tag  := line[1];
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
        if   tag = 'C' then
            Add(ops, rec(type := "const", id := nums[1]));
        elif tag = 'I' then
            Add(ops, rec(type := "inv",   a  := nums[1]));
        elif tag = 'M' then
            Add(ops, rec(type := "mul",   a  := nums[1], b := nums[2]));
        fi;
    od;
    CloseStream(file);
    return ops;
end;

# 3. File path
STRUCTURE_FILE := "structure.txt";
ops := ReadStructure(STRUCTURE_FILE);

# 4. Build leaf values
# const IDs in structure.txt are P-indices (C++ ordering)
# pToG converts them to G-indices so SLP runs in pc-group
maxLeafId  := Maximum(List(Filtered(ops, o -> o.type = "const"), o -> o.id));
leafValues := List([0 .. maxLeafId], i -> els[pToG[i + 1]]);

# 5. Evaluate SLP
EvaluateSLP := function(ops, leafValues)
    local values, i, op, inValA, inValB;
    values := ListWithIdenticalEntries(Length(ops), fail);
    for i in [1 .. Length(ops)] do
        op := ops[i];
        if op.type = "const" then
            values[i] := leafValues[op.id + 1];
        elif op.type = "inv" then
            values[i] := values[op.a + 1] ^ -1;
        else
            inValA    := values[op.a + 1];
            inValB    := values[op.b + 1];
            values[i] := inValA * inValB;
        fi;
    od;
    return values[Length(values)];
end;

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
for i in [1 .. warmupReps] do
    EvaluateSLP(ops, leafValues);
od;

GASMAN("collect");
start := Runtime();
for i in [1 .. repetitions] do
    result := EvaluateSLP(ops, leafValues);
od;
finish := Runtime();

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
while Length(fracStr) > 0 and fracStr[Length(fracStr)] = '0' do
    fracStr := fracStr{[1..Length(fracStr)-1]};
od;

if Length(fracStr) = 0 then
    timeStr := String(intPart);
else
    timeStr := Concatenation(String(intPart), ".", fracStr);
fi;

Print("\n--- GAP Results ---\n");
Print("Result index in P (0-based): ", resultIndexInP, "\n");
Print("Total time (ms): ", totalTime, "\n");
Print("Time per repetition (ms): ", timePerRep, "\n");

# 8. Write CSV


filepath := "benchmark_results.csv";

safeGroup := ReplacedString(GROUP_EXPR, ",", "_");
safeGroup := ReplacedString(safeGroup, "(", "");
safeGroup := ReplacedString(safeGroup, ")", "");

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
