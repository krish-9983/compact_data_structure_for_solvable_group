CayleyTableFor := function(H, G)
    local els, n, prodIdx, i, j, Sep, humanStr, parseStr, globalEls, els_global0, prodGlobal;

    els := Elements(H);
    n := Length(els);
    prodIdx := MultiplicationTable(H);

    globalEls := Elements(G);
    els_global0 := List(els, e -> Position(globalEls, e) - 1);

    prodGlobal := List([1..n], i -> []);

    for i in [1..n] do  
        for j in [1..n] do
            Add(prodGlobal[i], els_global0[prodIdx[i][j]]);

    
od;
od;

    Sep := function(j,n)
        if j = n then
            return "\n";
        else
            return "\t";
        fi;
    end;

    humanStr := Concatenation("\nCayley table for subgroup of order ", String(n), ":\n    ");
    for j in [1..n] do
        humanStr := Concatenation(humanStr, String(els_global0[j]), Sep(j,n));
    od;

    for i in [1..n] do
        humanStr := Concatenation(humanStr, String(els_global0[i]), "\t");
        for j in [1..n] do
            humanStr := Concatenation(humanStr, String(prodGlobal[i][j]), Sep(j,n));
        od;
    od;
    humanStr := Concatenation(humanStr, "\n");


    parseStr := Concatenation(String(n), " ");
    for i in [1..n] do
        for j in [1..n] do
            if i = n and j = n then
                parseStr := Concatenation(parseStr, String(prodGlobal[i][j]), "\n");
            else
                parseStr := Concatenation(parseStr, String(prodGlobal[i][j]), " ");
            fi;
        od;
    od;

    return rec(human := humanStr, parse := parseStr,
               elements := els_global0, table := prodGlobal);
end;


appendTextToFile := function(filename, text)
    local f;
    f := OutputTextFile(filename, true); 
    WriteAll(f, text);
    CloseStream(f);
end;


clearFile := function(filename)
    local f;
    f := OutputTextFile(filename, false);
    CloseStream(f);
end;




readablePath := "C:/Users/krish/Downloads/advanced_algorithm/goupCompactDS/cayley_readable.txt";
parsePath := "C:/Users/krish/Downloads/advanced_algorithm/goupCompactDS/cayley_parseable.txt";

clearFile(readablePath);
clearFile(parsePath);

G := SmallGroup(64,1);
cs := CompositionSeries(G);
out := rec();
for idx in [1..Length(cs)] do
    Print("=== Subgroup G", idx, " ( order ", Size(cs[idx]), ") === \n");
    out := CayleyTableFor(cs[idx], G);
    appendTextToFile(readablePath, out.human);
    appendTextToFile(parsePath, out.parse);
od;

Print("Appended all subgroups to:\n  readable: ", readablePath, "\n  parseable: ", parsePath, "\n");



#Read("C:/Users/krish/Downloads/advanced_algorithm/goupCompactDS/extractGroupInfo.g");