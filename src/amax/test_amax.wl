(* Test generalized a-maximization *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "FindCharges.wl"}]];

(* Test 1: SU(2) + 4 fundamentals (Nf=2, known result) *)
Print["=== Test 1: SU(2) + 4 fund ==="];
test1 = FindCharges[
  3, 2,  (* dimG=3, Tadj=2 *)
  {{"q", 1/2, 2, 4}},
  {}
];
Print[test1];

(* Test 2: SU(2) + 1 adj + 2 fund (CMNS seed 1) *)
Print["\n=== Test 2: SU(2) + 1 adj + 2 fund ==="];
test2 = FindCharges[
  3, 2,
  {{"q", 1/2, 2, 2}, {"phi", 2, 3, 1}},
  {}
];
Print[test2];

(* Test 3: SU(3) + 1 adj + 1 fund + 1 antifund (CMNS seed 3) *)
Print["\n=== Test 3: SU(3) + 1 adj + 1f + 1fb ==="];
test3 = FindCharges[
  8, 3,
  {{"q", 1/2, 3, 1}, {"qb", 1/2, 3, 1}, {"phi", 3, 8, 1}},
  {}
];
Print[test3];

(* Test 4: G2 + 1 adj + 1 fund (CMNS seed 10) *)
Print["\n=== Test 4: G2 + 1 adj + 1 fund ==="];
test4 = FindCharges[
  14, 4,
  {{"q", 1, 7, 1}, {"phi", 4, 14, 1}},
  {}
];
Print[test4];
