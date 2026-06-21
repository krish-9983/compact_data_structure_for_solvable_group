groups := [
  "SmallGroup(64,267)",
  "SmallGroup(72,41)",
  "SmallGroup(80,49)",
  "SmallGroup(96,227)",
  "SmallGroup(108,15)",
  "SmallGroup(112,41)",
  "SmallGroup(125,3)",
  "SmallGroup(128,928)",
  "SmallGroup(144,182)",
  "SmallGroup(162,10)",
  "SmallGroup(168,43)",
  "SmallGroup(192,1493)",
  "SmallGroup(200,40)",
  "SmallGroup(216,87)",
  "SmallGroup(224,63)",
  "SmallGroup(243,25)",
  "SmallGroup(252,26)",
  "SmallGroup(256,26675)",
  "SmallGroup(270,22)",
  "SmallGroup(288,1025)",
  "SmallGroup(300,25)",
  "SmallGroup(320,1635)",
  "SmallGroup(324,160)",
  "SmallGroup(336,118)",
  "SmallGroup(343,3)",
  "SmallGroup(378,6)",
  "SmallGroup(384,5765)",
  "SmallGroup(392,28)",
  "SmallGroup(400,206)",
  "SmallGroup(432,520)",
  "SmallGroup(448,179)",
  "SmallGroup(486,61)",
  "SmallGroup(500,27)",
  "SmallGroup(512,10494)",
  "SmallGroup(540,59)",
  "SmallGroup(576,8654)",
  "SmallGroup(600,55)",
  "SmallGroup(625,14)",
  "SmallGroup(648,259)",
  "SmallGroup(672,172)",
  "SmallGroup(686,5)",
  "SmallGroup(729,34)",
  "SmallGroup(750,15)",
  "SmallGroup(756,20)",
  "SmallGroup(768,1083755)",
  "SmallGroup(800,116)",
  "SmallGroup(810,65)",
  "SmallGroup(864,700)",
  "SmallGroup(882,7)",
  "SmallGroup(896,211)",
  "SmallGroup(900,95)",
  "SmallGroup(972,66)",
  "SmallGroup(1000,86)",
  "SmallGroup(1029,12)",
  "SmallGroup(1080,104)",
  "SmallGroup(1125,10)",
  "SmallGroup(1152,157)",
  "SmallGroup(1176,41)",
  "SmallGroup(1215,15)",
  "SmallGroup(1250,6)",
  "SmallGroup(1296,536)",
  "SmallGroup(1323,4)",
  "SmallGroup(1372,35)",
  "SmallGroup(1440,523)",
  "SmallGroup(1458,30)",
  "SmallGroup(1500,28)",
  "SmallGroup(1536,408544)",
  "SmallGroup(1568,41)",
  "SmallGroup(1620,63)",
  "SmallGroup(1680,954)",
  "SmallGroup(1701,28)",
  "SmallGroup(1728,928)",
  "SmallGroup(1750,10)",
  "SmallGroup(1944,368)",
  "SmallGroup(2000,654)",
  "SmallGroup(2187,289)",
  "SmallGroup(3125,5)",
  "SmallGroup(3645,10)",
  "DirectProduct(SmallGroup(1029,12),CyclicGroup(2))",
  "DirectProduct(SmallGroup(1080,104),CyclicGroup(2))",
  "DirectProduct(SmallGroup(1125,10),CyclicGroup(2))",
  "DirectProduct(SmallGroup(1152,157),CyclicGroup(2))",
  "DirectProduct(SmallGroup(1176,41),CyclicGroup(2))",
  "DirectProduct(SmallGroup(1215,15),CyclicGroup(2))",
  "DirectProduct(SmallGroup(1250,6),CyclicGroup(2))",
  "DirectProduct(SmallGroup(1296,536),CyclicGroup(2))",
  "DirectProduct(SmallGroup(1323,4),CyclicGroup(2))",
  "DirectProduct(SmallGroup(1372,35),CyclicGroup(2))",
  "DirectProduct(SmallGroup(1440,523),CyclicGroup(2))",
  "DirectProduct(SmallGroup(1458,30),CyclicGroup(2))",
  "DirectProduct(SmallGroup(1500,28),CyclicGroup(2))",
  "DirectProduct(SmallGroup(1536,408544),CyclicGroup(2))",
  "DirectProduct(SmallGroup(1568,41),CyclicGroup(2))",
  "DirectProduct(SmallGroup(1620,63),CyclicGroup(2))",
  "DirectProduct(SmallGroup(1680,954),CyclicGroup(2))",
  "DirectProduct(SmallGroup(1701,28),CyclicGroup(2))",
  "DirectProduct(SmallGroup(1728,928),CyclicGroup(2))",
  "DirectProduct(SmallGroup(1750,10),CyclicGroup(2))",
  "DirectProduct(SmallGroup(1944,368),CyclicGroup(2))",
  "DirectProduct(SmallGroup(2000,654),CyclicGroup(2))",
  "DirectProduct(SmallGroup(2187,289),CyclicGroup(2))"
];

ok_count   := 0;
fail_count := 0;
ns_count   := 0;

for expr in groups do
    res := CALL_WITH_CATCH( EvalString, [ expr ] );
    if res[1] = false then
        Print( "FAIL         ", expr, "\n" );
        fail_count := fail_count + 1;
    elif not IsGroup( res[2] ) then
        Print( "FAIL         ", expr, "\n" );
        fail_count := fail_count + 1;
    elif not IsSolvable( res[2] ) then
        Print( "NOT_SOLVABLE ", expr, "\n" );
        ns_count := ns_count + 1;
    else
        G := res[2];
        Print( "OK  order=", Size(G),
               "  ", expr, "\n" );
        ok_count := ok_count + 1;
    fi;
od;

Print( "\n============================================================\n" );
Print( "OK: ", ok_count,
       "  FAIL: ", fail_count,
       "  NOT_SOLVABLE: ", ns_count, "\n" );
Print( "============================================================\n" );
QUIT_GAP(0);
