(* Test OperatorSpectrum *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "OperatorSpectrum.wl"}]];

(* Test 1: SU(2) + 8 fund, W=0, R[q]=0.5 *)
Print["=== Test 1: SU(2) 8x fund, W=0 ==="];
result1 = EnumerateOperators[
  {{"r0", 8, "pseudo-real", 1/2}},
  {}
];
Print["n_all=", result1["n_all"],
  " n_relevant=", result1["n_relevant"],
  " n_superrel=", result1["n_superrelevant"]];
Print["Relevant ops (first 5):"];
Do[Print["  ", r], {r, Take[result1["relevant_with_R"], Min[5, Length[result1["relevant_with_R"]]]]}];

(* Test 2: SU(2) + 8 fund, W=M1*r01*r02, R[q1]=R[q2]=0.54, R[q3-8]=0.487, R[M]=0.919 *)
Print["\n=== Test 2: SU(2) 8x fund + 1 singlet, W=M1*r01*r02 ==="];
result2 = EnumerateOperators[
  {{"r0", 8, "pseudo-real", 1/2},
   {"M", 1, "singlet", 0.919}},
  {"M1*r01*r02"}
];
Print["n_all=", result2["n_all"],
  " n_relevant=", result2["n_relevant"],
  " n_superrel=", result2["n_superrelevant"]];
Print["Relevant ops (first 10):"];
Do[Print["  ", r], {r, Take[result2["relevant_with_R"], Min[10, Length[result2["relevant_with_R"]]]]}];

(* Test 3: SU(2) + 2 fund + 1 adj, W=0 (CMNS seed 1) *)
Print["\n=== Test 3: SU(2) 2f + 1adj, W=0 ==="];
result3 = EnumerateOperators[
  {{"q", 2, "pseudo-real", 0.4767},
   {"phi", 1, "real", 0.2617}},
  {}
];
Print["n_all=", result3["n_all"],
  " n_relevant=", result3["n_relevant"],
  " n_superrel=", result3["n_superrelevant"]];
Print["Relevant ops:"];
Do[Print["  ", r], {r, result3["relevant_with_R"]}];
