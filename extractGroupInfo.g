# extractGroupInfo.g  —  Stage 1: GAP file for composition series analysis and cayley table analysis

# HOW TO CALL
#   gap -q -b -c 'GROUP_EXPR:="SmallGroup(96,227)";' extractGroupInfo.g

# FILES WRITTEN
#   cayley_readable.txt   — human-readable tables, useful for debugging
#   cayley_parseable.txt  — machine-readable format read by preparedsa.cpp
#   group_order.txt       — single integer |G|, used by run_all.sh to bound
#                           number of random elements for structure_generator.c file

#   I convert the polycyclic group to a permutation group via
#   IsomorphismPermGroup before doing anything. These permutation groups have a canonical
#   element ordering that GAP reproduces identically across all sessions. 

# THE GLOBAL INDEX SCHEME 
#   Every entry in every subgroup's Cayley table is stored as a 0-based index
#   into G's full element list, NOT as a local H-index. This choice
#   means preparedsa and querycds never need a local-to-global translation map
#   anywhere.
# =============================================================================


#This function computes the Cayley table of subgroup H ( in global G-indices)
#els_global0: it translates each H-element from its local position in H to its 0-based position in G.

CayleyTableFor := function(H, G, readablePath, parsePath)
    local els, n, prodIdx, i, j, val,
          globalEls, els_global0,
          fHuman, fParse;

    #all elements of H, where each element is a permutation
    els  := Elements(H);
    n    := Length(els);

    # MultiplicationTable returns an n x n matrix of LOCAL 1-based indices.
    # We need to translate those to global 0-based G-indices below.
    prodIdx   := MultiplicationTable(H);
    globalEls := Elements(G);


    #It converts from GAP's local 1-based to our global0-based ordering
    els_global0 := List(els,
        e -> PositionCanonical(globalEls, e) - 1);

   #file opening
    fHuman := OutputTextFile(readablePath, true);
    fParse := OutputTextFile(parsePath, true);

    # human-readable header
    WriteAll(fHuman, Concatenation(
        "\nCayley table for subgroup of order ", String(n), ":\n    "));

    for j in [1..n] do
        WriteAll(fHuman, String(els_global0[j]));
        if j = n then WriteAll(fHuman, "\n");
        else          WriteAll(fHuman, "\t");
        fi;
    od;

    # ---- parseable format: "n v00 v01 ... v(n-1)(n-1)" ---------------
    # preparedsa reads the order n first, then n*n values. The whole table
    # is on one line, space-separated, with a newline at the very end.
    WriteAll(fParse, String(n));
    WriteAll(fParse, " ");

    for i in [1..n] do
        WriteAll(fHuman, String(els_global0[i]));
        WriteAll(fHuman, "\t");

        for j in [1..n] do
            # prodIdx[i][j] is a local 1-based H-index.
            # els_global0[ .... ] converts it to the global 0-based G-index.
            val := els_global0[ prodIdx[i][j] ];

            WriteAll(fHuman, String(val));
            if j = n then WriteAll(fHuman, "\n");
            else          WriteAll(fHuman, "\t");
            fi;

            WriteAll(fParse, String(val));
            # Newline only after the very last value of this subgroup's table.
            if i = n and j = n then WriteAll(fParse, "\n");
            else                     WriteAll(fParse, " ");
            fi;
        od;
    od;

    CloseStream(fHuman);
    CloseStream(fParse);
end;


# This clears the file before writing.
clearFile := function(filename)
    local f;
    f := OutputTextFile(filename, false);   # false = truncate, not append
    CloseStream(f);
end;


# ---- Main Part ----

readablePath := "cayley_readable.txt";
parsePath    := "cayley_parseable.txt";

# IT always clear before writing so data from previous group run clear and don't produce wrong data.
clearFile(readablePath);
clearFile(parsePath);

# GROUP_EXPR is set by run_all.sh before calling this script. if not set then produce error.
if not IsBound(GROUP_EXPR) then
    Error("Use: gap -q -c 'GROUP_EXPR:=\"SmallGroup(96,227)\";' extractGroupInfo.g");
fi;

P := EvalString(GROUP_EXPR);

Print("Group loaded: ", GROUP_EXPR, "\n");
Print("Order: ", Size(P), "\n");

#writing order of Group
fOrder := OutputTextFile("group_order.txt", false);
WriteAll(fOrder, String(Size(P)));
WriteAll(fOrder, "\n");
CloseStream(fOrder);

#following Convert to permutation group for canonical, reproducible element ordering.
# This also makes Elements(G) and MultiplicationTable(H) deterministic across GAP sessions.

iso := IsomorphismPermGroup(P);
G   := Image(iso);
cs  := CompositionSeries(G);

# Write Cayley tables for ALL subgroups in the composition series.
for idx in [1..Length(cs)] do
    Print("=== Subgroup G", idx, " ( order ", Size(cs[idx]), ") ===\n");
    CayleyTableFor(cs[idx], G, readablePath, parsePath);
od;

Print("Appended all subgroups to:\n");
Print("  readable : ", readablePath, "\n");
Print("  parseable: ", parsePath,    "\n");

Print("Ended", "\n");

