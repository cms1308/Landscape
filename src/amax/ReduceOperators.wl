(* ::Package:: *)
(*
  ReduceOperators: Apply F-term reduction to a list of operator monomials.

  Takes externally-generated monomial list (from Python orbit-aware PE pipeline)
  and applies GroebnerBasis reduction from superpotential F-terms.

  ReduceAndClassify[monomialStrings, fieldNames, w, rCharges, maxR]

  Args:
    monomialStrings: list of monomial strings {"r01*r03", "M01*r03*r04", ...}
    fieldNames: list of all field name strings {"r01", "r02", ..., "M01"}
    w: list of superpotential term strings {"M01*r01*r02"}
    rCharges: Association mapping field name -> R-charge
    maxR: keep operators with R < maxR (default 2)

  Returns: Association with
    "reduced": list of surviving monomial strings
    "relevant": monomials with R < 2
    "superrelevant": monomials with R < 4/3
    "flippable": monomials with 2/3 < R < 4/3
    "unitarity_violating": monomials with R <= 2/3
    "relevant_with_R": {monomial, R} pairs
*)

ReduceAndClassify[monomialStrings_List, fieldNames_List, w_List, rCharges_Association, maxR_:2] :=
 Module[{
   flatFields, w2, fterms, ringRelations,
   monomials, reduced, opcharge,
   relevant, superrelevant, flippable, unitarityViol},

  (* Convert field names to symbols *)
  flatFields = ToExpression /@ fieldNames;

  (* R-charge function *)
  Do[opcharge[ToExpression[k]] = rCharges[k], {k, Keys[rCharges]}];
  opcharge[expr_Times] := Plus @@ (opcharge /@ (List @@ expr));
  opcharge[Power[base_, n_Integer]] := n * opcharge[base];

  (* Parse superpotential *)
  w2 = If[w === {} || w === {""}, {},
    ToExpression /@ Select[w, # =!= "" &]];

  (* Convert monomial strings to expressions *)
  monomials = ToExpression /@ monomialStrings;

  (* F-term reduction *)
  If[Length[w2] > 0,
    fterms = Table[D[Total[w2], v], {v, flatFields}];
    fterms = DeleteCases[fterms, 0];
    If[Length[fterms] > 0,
      ringRelations = GroebnerBasis[fterms, flatFields];
      (* Reduce each monomial *)
      reduced = DeleteDuplicates[DeleteCases[
        (#[[2]] /. Times[c_?NumberQ, rest__] :> Times[rest]) & /@
          PolynomialReduce[monomials, ringRelations, flatFields],
        0
      ]];
      (* Keep only monomials (not sums — those arise from incomplete reduction) *)
      reduced = Select[reduced, Head[#] =!= Plus &];
      ,
      reduced = monomials;
    ],
    reduced = monomials;
  ];

  (* Remove duplicates and filter by R *)
  reduced = DeleteDuplicates[reduced];
  reduced = Select[reduced, opcharge[#] < maxR &];

  (* Classify *)
  relevant = Select[reduced, opcharge[#] < 2 &];
  superrelevant = Select[reduced, opcharge[#] < 4/3 &];
  flippable = Select[reduced, 2/3 + 10^-10 < opcharge[#] < 4/3 &];
  unitarityViol = Select[reduced, opcharge[#] <= 2/3 + 10^-10 &];

  <|
    "reduced" -> (ToString[InputForm[#]] & /@ reduced),
    "relevant" -> (ToString[InputForm[#]] & /@ relevant),
    "superrelevant" -> (ToString[InputForm[#]] & /@ superrelevant),
    "flippable" -> (ToString[InputForm[#]] & /@ flippable),
    "unitarity_violating" -> (ToString[InputForm[#]] & /@ unitarityViol),
    "relevant_with_R" -> ({ToString[InputForm[#]], N[opcharge[#], 10]} & /@ relevant),
    "n_reduced" -> Length[reduced],
    "n_relevant" -> Length[relevant],
    "n_superrelevant" -> Length[superrelevant],
    "n_flippable" -> Length[flippable],
    "n_unitarity_violating" -> Length[unitarityViol]
  |>
]
