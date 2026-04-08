(* ::Package:: *)
(*
  OperatorSpectrum: enumerate GIOs, apply F-term reduction, classify.

  EnumerateOperators[repInfo, nCopies, w, rcharges]

  Constructs all gauge-invariant monomials from:
  1. Gauge-invariant pairs (mesons) from PE degree-2 data
     - Real reps: symmetric pairing (i <= j)
     - Pseudo-real reps: antisymmetric pairing (i < j)
     - Mixed pairs: all combinations
  2. Singlet fields (M, X) and their products with gauge invariants
  3. Products of mesons (quartic, sextic, ...) up to R < 2
  4. Apply F-term reduction via GroebnerBasis
  5. Classify: relevant (R<2), super-relevant (R<4/3)
*)

EnumerateOperators[fieldInfo_List, w_List, maxR_:2] :=
 Module[{
   fields, flatFields, singletFields, gaugeFields,
   w2, fterms, ringRelations,
   mesons, quartic, allGIO, gio, reduced,
   opcharge, relevant, superrelevant, flippable},

  (* fieldInfo: list of {name, nCopies, reality, rChargeOrList}
     name: field name prefix (e.g. "r0", "M0")
     nCopies: number of copies
     reality: "real", "pseudo-real", "complex", or "singlet"
     rChargeOrList: single number (all copies same) OR list of per-copy R-charges *)

  (* Generate field variable names *)
  fields = Table[
    Table[ToExpression[fieldInfo[[i, 1]] <> ToString[j]],
      {j, fieldInfo[[i, 2]]}],
    {i, Length[fieldInfo]}
  ];
  flatFields = Flatten[fields];

  singletFields = Flatten[Table[
    If[fieldInfo[[i, 4]] == "singlet" || fieldInfo[[i, 3]] == "singlet",
      fields[[i]], {}],
    {i, Length[fieldInfo]}
  ]];

  gaugeFields = Flatten[Table[
    If[fieldInfo[[i, 3]] =!= "singlet",
      fields[[i]], {}],
    {i, Length[fieldInfo]}
  ]];

  (* R-charge function: assign to each field variable *)
  Do[
    With[{rc = fieldInfo[[i, 4]]},
      If[ListQ[rc],
        (* Per-copy R-charges *)
        Do[opcharge[fields[[i, j]]] = rc[[j]], {j, fieldInfo[[i, 2]]}],
        (* Single R-charge for all copies *)
        Do[opcharge[fields[[i, j]]] = rc, {j, fieldInfo[[i, 2]]}]
      ]
    ],
    {i, Length[fieldInfo]}
  ];
  (* For composite operators: expand into base variables and sum *)
  opcharge[expr_Times] := Plus @@ (opcharge /@ (List @@ expr));
  opcharge[Power[base_, n_Integer]] := n * opcharge[base];

  (* === 1. Construct gauge-invariant mesons === *)
  mesons = {};

  (* Pairs within same rep type *)
  Do[
    If[fieldInfo[[i, 3]] === "singlet", Continue[]];
    With[{fs = fields[[i]], reality = fieldInfo[[i, 3]]},
      If[reality === "pseudo-real",
        (* Antisymmetric pairing: i < j *)
        mesons = Join[mesons,
          (Times @@ #) & /@ Subsets[fs, {2}]],
        (* Real: symmetric pairing: i <= j *)
        mesons = Join[mesons,
          (Times @@ #) & /@ Join[
            Subsets[fs, {2}],
            {#^2} & /@ fs
          ]]
      ]
    ],
    {i, Length[fieldInfo]}
  ];

  (* Pairs between different gauge-charged rep types *)
  Do[
    If[fieldInfo[[i, 3]] === "singlet" || fieldInfo[[j, 3]] === "singlet",
      Continue[]];
    mesons = Join[mesons,
      Flatten[Outer[Times, fields[[i]], fields[[j]]]]
    ],
    {i, Length[fieldInfo]}, {j, i + 1, Length[fieldInfo]}
  ];

  (* === 2. Quartic: products of two mesons === *)
  quartic = DeleteDuplicates[
    Times @@@ Tuples[mesons, 2]
  ];

  (* === 3. Full GIO list: singlets + mesons + quartic + products === *)
  allGIO = Join[
    singletFields,                              (* singlet fields *)
    mesons,                                     (* mesons *)
    Flatten[Outer[Times, singletFields, mesons]],  (* singlet * meson *)
    Flatten[Outer[Times, singletFields, singletFields]], (* singlet * singlet *)
    quartic,                                    (* meson * meson *)
    Flatten[Outer[Times, singletFields, quartic]]  (* singlet * quartic *)
  ];

  (* Remove duplicates and zero *)
  allGIO = DeleteDuplicates[DeleteCases[allGIO, 0]];

  (* Filter by R-charge: only keep R < maxR *)
  allGIO = Select[allGIO, opcharge[#] < maxR &];

  (* === 4. F-term reduction === *)
  w2 = If[w === {} || w === {""}, {},
    ToExpression /@ Select[w, # =!= "" &]];

  If[Length[w2] > 0,
    fterms = Table[D[Total[w2], v], {v, flatFields}];
    fterms = DeleteCases[fterms, 0];
    If[Length[fterms] > 0,
      ringRelations = GroebnerBasis[fterms, flatFields];
      (* Reduce each GIO modulo ring relations *)
      gio = DeleteDuplicates[DeleteCases[
        (#[[2]] /. Times[c_?NumberQ, rest__] :> Times[rest]) & /@
          PolynomialReduce[allGIO, ringRelations, flatFields],
        0
      ]];
      (* Keep only monomials (not sums) *)
      gio = Select[gio, Head[#] =!= Plus &];
      ,
      gio = allGIO;
    ],
    gio = allGIO;
  ];

  (* === 5. Classify === *)
  relevant = Select[gio, opcharge[#] < 2 &];
  superrelevant = Select[gio, opcharge[#] < 4/3 &];
  flippable = Select[gio, 2/3 < opcharge[#] < 4/3 &];

  (* Return *)
  <|
    "all_gio" -> (ToString[InputForm[#]] & /@ gio),
    "relevant" -> (ToString[InputForm[#]] & /@ relevant),
    "superrelevant" -> (ToString[InputForm[#]] & /@ superrelevant),
    "flippable" -> (ToString[InputForm[#]] & /@ flippable),
    "relevant_with_R" -> ({ToString[InputForm[#]], N[opcharge[#], 10]} & /@ relevant),
    "n_all" -> Length[gio],
    "n_relevant" -> Length[relevant],
    "n_superrelevant" -> Length[superrelevant],
    "n_flippable" -> Length[flippable]
  |>
]
