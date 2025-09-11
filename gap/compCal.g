n := 4;
mt := [
    [1,2,3,4],
    [2,3,4,1],
    [3,4,1,2],
    [4,1,2,3]
];

perms := List([1..n], i -> PermList(mt[i]));
G := Group(perms);
if Size(G) <> n then Error("Regular action failed"); fi;

series := CompositionSeries(G);
m := Length(series) - 1;
Display(List(series, Size));

permIndex := rec();
for a in [1..n] do
  permIndex.(String(perms[a])) := a;
od;

IdsOfSubgroup := function(H)
  local ids, x;
  ids := [];
  for x in Elements(H) do
    Add(ids, permIndex.(String(x)));
  od;
  Sort(ids);
  return ids;
end;

chainIds := List(series, H -> IdsOfSubgroup(H));

layerInfo := [];
for i in [1..m] do
  Ki  := series[i];
  Kip := series[i+1];
  pi  := NaturalHomomorphismByNormalSubgroup(Ki, Kip);
  Qi  := Image(pi);
  ki  := Size(Qi);

  genQ := One(Qi);
  for x in GeneratorsOfGroup(Qi) do
    if Order(x) = ki then genQ := x; break; fi;
  od;
  if genQ = One(Qi) then
    for x in Elements(Qi) do if Order(x) = ki then genQ := x; break; fi; od;
  fi;
  if genQ = One(Qi) then Error("Could not find generator of quotient at layer ", i); fi;

  pre := PreImagesRepresentative(pi, genQ);

  Add(layerInfo, rec(
    i      := i,
    sizeKi := Size(Ki),
    sizeKip:= Size(Kip),
    k      := ki,
    KiIds  := IdsOfSubgroup(Ki),
    KipIds := IdsOfSubgroup(Kip),
    g0Id   := permIndex.(String(pre))
  ));
od;

# === absolute file paths (edit these to match your system!) ===
# === absolute file paths (edit these to match your system!) ===
compFile := "C:/Users/krish/Downloads/advanced algorithm/thesis/gap/comp_series.txt";
layerFile := "C:/Users/krish/Downloads/advanced algorithm/thesis/gap/layers.csv";

PrintTo(compFile, "# each line = 1-based IDs of Ki (K0 down to Km)\n");
for ids in chainIds do
  AppendTo(compFile, JoinStringsWithSeparator(List(ids, String), " "), "\n");
od;

PrintTo(layerFile,
  "# i,sizeKi,sizeKip,k,g0Id,|SEP|,KiIds...,|SEP|,KipIds...,|SEP|,reps...\n");

for reci in layerInfo do
  # meta fields first
  AppendTo(layerFile,
    reci.i, ",", reci.sizeKi, ",", reci.sizeKip, ",", reci.k, ",", reci.g0Id, ",|SEP|,");

  # Ki
  AppendTo(layerFile, JoinStringsWithSeparator(List(reci.KiIds, String), " "), ",|SEP|,");

  # Kip
  AppendTo(layerFile, JoinStringsWithSeparator(List(reci.KipIds, String), " "), ",|SEP|,");

  # reps (currently not computed, so leave blank)
  AppendTo(layerFile, "\n");
od;

