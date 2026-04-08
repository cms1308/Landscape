(* Run SU(2) 8x fund through depth 2 with full operator spectrum *)

Get[FileNameJoin[{Directory[], "src/amax/FindCharges.wl"}]];
Get[FileNameJoin[{Directory[], "src/amax/OperatorSpectrum.wl"}]];

(* === Helper: run one theory === *)
AnalyzeTheory[dimG_, Tadj_, repInfo_, w_, tag_] := Module[
  {result, fieldInfo, ops, rcharges},

  (* a-maximization *)
  result = FindCharges[dimG, Tadj, repInfo, w];

  If[result["consistency"] =!= "consistent",
    Print["  ", tag, ": ", result["consistency"]];
    Return[<|"tag" -> tag, "consistency" -> result["consistency"]|>]
  ];

  (* Build fieldInfo for OperatorSpectrum *)
  fieldInfo = Table[
    With[{name = repInfo[[i, 1]], nCop = repInfo[[i, 4]],
          T = repInfo[[i, 2]], dim = repInfo[[i, 3]]},
      {name, nCop,
       If[T === 0, "singlet",
         If[dim === 1, "singlet", "pseudo-real"]], (* SU(2): all non-singlet reps are pseudo-real *)
       N[result[name <> "1"], 15]}
    ],
    {i, Length[repInfo]}
  ];

  (* Enumerate operators *)
  ops = EnumerateOperators[fieldInfo, w];

  (* Check unitarity *)
  unitarityViol = Select[ops["relevant_with_R"],
    #[[2]] <= 2/3 + 10^-10 &];

  If[Length[unitarityViol] > 0,
    Print["  ", tag, ": operator decoupled (", Length[unitarityViol], " ops with R<=2/3)"];
    Return[<|"tag" -> tag, "consistency" -> "operator decoupled",
      "a" -> N[result["a"], 15], "c" -> N[result["c"], 15]|>]
  ];

  Print["  ", tag, ": consistent, a=", N[result["a"], 8],
    ", relevant=", ops["n_relevant"],
    ", superrel=", ops["n_superrelevant"],
    ", flippable=", ops["n_flippable"]];

  <|"tag" -> tag,
    "consistency" -> "consistent",
    "a" -> N[result["a"], 15],
    "c" -> N[result["c"], 15],
    "repInfo" -> repInfo,
    "w" -> w,
    "relevant" -> ops["relevant"],
    "superrelevant" -> ops["superrelevant"],
    "flippable" -> ops["flippable"],
    "relevant_with_R" -> ops["relevant_with_R"],
    "fieldInfo" -> fieldInfo,
    "result" -> result
  |>
];

(* === Helper: generate deformations from a theory === *)
GenerateDeformations[theory_, parentIdx_] := Module[
  {deformations = {}, relevant, flippable, repInfo, w, op, newW, newRepInfo, singletName},

  If[theory["consistency"] =!= "consistent", Return[{}]];

  relevant = theory["relevant"];
  flippable = theory["flippable"];
  repInfo = theory["repInfo"];
  w = theory["w"];

  (* Direct deformations *)
  Do[
    newW = Append[w, op];
    AppendTo[deformations, <|
      "repInfo" -> repInfo,
      "w" -> newW,
      "type" -> "direct",
      "op" -> op,
      "parent" -> parentIdx
    |>],
    {op, relevant}
  ];

  (* Flip deformations *)
  Do[
    singletName = "M" <> ToString[Length[repInfo]];
    newRepInfo = Append[repInfo, {singletName, 0, 1, 1}];
    newW = Append[w, singletName <> "1*" <> op];
    AppendTo[deformations, <|
      "repInfo" -> newRepInfo,
      "w" -> newW,
      "type" -> "flip",
      "op" -> op,
      "parent" -> parentIdx
    |>],
    {op, flippable}
  ];

  deformations
];

(* === Main: SU(2) 8x fund === *)
Print["========================================"];
Print["  SU(2) + 8 fund, depth 0 to 2"];
Print["========================================"];

dimG = 3; Tadj = 2;
baseRepInfo = {{"r0", 1/2, 2, 8}};

(* --- Depth 0 --- *)
Print["\n--- Depth 0 ---"];
depth0 = {AnalyzeTheory[dimG, Tadj, baseRepInfo, {}, "d0"]};
Print["Depth 0: ", Length[Select[depth0, #["consistency"] === "consistent" &]], " consistent"];

(* --- Depth 1 --- *)
Print["\n--- Depth 0 -> 1 ---"];
cands1 = Flatten[Table[
  GenerateDeformations[t, i],
  {i, Length[depth0]}, {t, {depth0[[i]]}}
], 2];
Print["  ", Length[cands1], " candidates before dedup"];

(* Deduplicate: for W=0 seed, all mesons equivalent, so keep one direct + one flip *)
(* Simple dedup: group by (type, sorted W structure) *)
dedupKeys1 = DeleteDuplicates[
  Table[{cands1[[i, "type"]],
    Sort[StringSplit[#, "*"] & /@ cands1[[i, "w"]]]}, {i, Length[cands1]}]
];
Print["  ", Length[dedupKeys1], " after dedup"];

(* Run unique candidates *)
depth1 = {};
uniqueIdx = DeleteDuplicatesBy[Range[Length[cands1]],
  {cands1[[#, "type"]],
   Sort[Length /@ (StringSplit[#, "*"] & /@ cands1[[#, "w"]])]} &
];
Do[
  c = cands1[[idx]];
  tag = "d1_" <> c["type"] <> "_" <> ToString[idx];
  t = AnalyzeTheory[dimG, Tadj, c["repInfo"], c["w"], tag];
  If[t["consistency"] === "consistent",
    t = Append[t, "deformType" -> c["type"]];
    t = Append[t, "deformOp" -> c["op"]];
    AppendTo[depth1, t]
  ],
  {idx, uniqueIdx}
];
Print["Depth 1: ", Length[depth1], " consistent"];

(* --- Depth 2 --- *)
Print["\n--- Depth 1 -> 2 ---"];
cands2 = Flatten[Table[
  GenerateDeformations[depth1[[i]], i],
  {i, Length[depth1]}
]];
Print["  ", Length[cands2], " candidates before dedup"];

(* Dedup for depth 2 *)
uniqueIdx2 = DeleteDuplicatesBy[Range[Length[cands2]],
  Sort[StringSplit[#, "*"] & /@ cands2[[#, "w"]]] &
];
Print["  ", Length[uniqueIdx2], " after dedup"];

depth2 = {};
Do[
  c = cands2[[idx]];
  tag = "d2_" <> c["type"] <> "_" <> ToString[idx];
  t = AnalyzeTheory[dimG, Tadj, c["repInfo"], c["w"], tag];
  If[t["consistency"] === "consistent",
    AppendTo[depth2, t]
  ],
  {idx, uniqueIdx2}
];
Print["Depth 2: ", Length[depth2], " consistent"];

(* Summary *)
Print["\n========================================"];
Print["  Summary"];
Print["========================================"];
Print["Depth 0: ", Length[Select[depth0, #["consistency"] === "consistent" &]]];
Print["Depth 1: ", Length[depth1]];
Print["Depth 2: ", Length[depth2]];
Print["Total: ", Length[Select[depth0, #["consistency"] === "consistent" &]] + Length[depth1] + Length[depth2]];
