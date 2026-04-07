(* ::Package:: *)
(*
  Generalized a-maximization for arbitrary simple gauge groups.

  FindCharges[dimG, Tadj, repInfo, w]

  Arguments:
    dimG    : dimension of gauge group (integer)
    Tadj    : T(adjoint) = h^v (rational)
    repInfo : list of {name, T_R, dimR, nCopies} for each matter rep type
              Singlets: {"X", 0, 1, nX}, {"M", 0, 1, nM}
              Gauge-charged: {"q", 1/2, Nc, nq}, etc.
    w       : {} or {"term1", "term2", ...} superpotential monomials
              Fields named as name<>index: q1, q2, M1, phi1, etc.

  Returns: Association with keys:
    "consistency", "method", "a", "c", per-field R-charges,
    "rational", "computation_time"
*)

FindCharges[dimG_Integer, Tadj_, repInfo_List, w_List] :=
 Block[{$MaxExtraPrecision = 100},
  Module[{
    q, names, mu, dims, nCopies, mat, flatMat, w2,
    abjSum, abjVar, abjVarBase,
    getMonomialCharge, termSum, v, baseVars, selectedVarBase, val,
    RRR, R, a, c, xVars, res, Hess, localmax,
    indepVars, bestSol, rcharge, central, rational, consistency,
    result, method, time, finalResult,
    singletNames, fterms, ringRelations},

   {time, finalResult} = AbsoluteTiming[

    (* Unpack repInfo *)
    names = repInfo[[All, 1]];
    mu = repInfo[[All, 2]];
    dims = repInfo[[All, 3]];
    nCopies = repInfo[[All, 4]];

    (* Singlet names: those with T = 0 *)
    singletNames = Select[names, mu[[Position[names, #][[1, 1]]]] == 0 &];

    w2 = If[w === {}, {}, ToExpression /@ w];
    method = <|"method" -> "Exactly solved"|>;

    (* Field variables *)
    mat = Table[
      ToExpression[names[[i]] <> ToString[j]],
      {i, Length[nCopies]}, {j, nCopies[[i]]}
    ];
    flatMat = Flatten[mat];

    Catch[

     (* ===== 1. ABJ Anomaly: Tadj + sum (q[field]-1)*T(rep) = 0 ===== *)
     abjSum = Tadj + Sum[
        (q[mat[[i, j]]] - 1) mu[[i]],
        {i, Length[nCopies]}, {j, nCopies[[i]]}
       ];

     If[Variables[abjSum] =!= {},
      abjVar = First[Variables[abjSum]];
      abjVarBase = abjVar /. q[sym_] :> sym;
      q[abjVarBase] = abjVar /. First[Solve[abjSum == 0, abjVar]];
     ];

     (* ===== 2. Superpotential: R[term] = 2 for each term ===== *)
     getMonomialCharge[monomial_] :=
      Total[Cases[Variables[monomial],
        var_ :> Exponent[monomial, var]*q[var]]];

     Do[
      termSum = Expand[getMonomialCharge[term]];
      If[termSum =!= 2,
       v = Variables[termSum];
       If[Length[v] == 0,
        Throw[<|"n" -> nCopies, "w" -> w,
          "consistency" -> "Superpotential not preserving R-symmetry"|>]
       ];
       baseVars = v /. q[sym_] :> sym;
       (* Prioritize singlet fields *)
       selectedVarBase = With[{sn = singletNames},
        If[AnyTrue[baseVars,
          Function[bv, AnyTrue[sn, StringStartsQ[ToString[bv], #] &]]],
         First[Select[baseVars,
           Function[bv, AnyTrue[sn, StringStartsQ[ToString[bv], #] &]]]],
         First[baseVars]
        ]];
       val = Solve[termSum == 2, q[selectedVarBase]];
       If[Length[val] == 0,
        Throw[<|"n" -> nCopies, "w" -> w,
          "consistency" -> "Superpotential not preserving R-symmetry"|>],
        q[selectedVarBase] = q[selectedVarBase] /. First[val];
       ];
      ],
      {term, w2}
     ];

     (* ===== 3. Central Charges ===== *)
     (* Tr R^3 = dimG + sum_i nCopies_i * (q[field]-1)^3 * dim(R_i) *)
     RRR = dimG + Sum[
        (q[mat[[i, j]]] - 1)^3 dims[[i]],
        {i, Length[nCopies]}, {j, nCopies[[i]]}
       ];
     R = dimG + Sum[
        (q[mat[[i, j]]] - 1) dims[[i]],
        {i, Length[nCopies]}, {j, nCopies[[i]]}
       ];
     a = Expand[3 (3 RRR - R)/32];
     c = Expand[(9 RRR - 5 R)/32];

     (* ===== 4. a-Maximization ===== *)
     xVars = Variables[a];

     If[Length[xVars] === 0,
      (* Fully determined *)
      Null,

      (* Check variable count *)
      If[Total[nCopies] - Length[w] - 1 =!= Length[xVars],
       (* Mismatch - may have redundant constraints or symmetry *)
       Null (* proceed anyway *)
      ];

      indepVars = Variables[xVars];

      If[Length[indepVars] > 0,
       (* Try exact solve *)
       res = Quiet[TimeConstrained[
         Simplify[TimeConstrained[
           Solve[Thread[D[a, {indepVars}] == 0], indepVars],
           600, $NSolve]],
         600, $NSolve]];

       If[res =!= $NSolve,
        res = Select[res, AllTrue[Values[#], NumericQ] &];
       ];

       (* Fallback: numerical *)
       If[res === $NSolve || Length[res] == 0,
        method = <|"method" -> "Numerically solved"|>;
        res = TimeConstrained[
          Quiet[Check[
            NSolve[Thread[D[a, {indepVars}] == 0], indepVars,
             WorkingPrecision -> 50],
            $Failed, {NSolve::infsolns}]],
          600, $Failed];
       ];

       If[res === $Failed,
        Throw[<|"n" -> nCopies, "w" -> w,
          "consistency" -> "a-maximization timeout"|>]
       ];

       (* Select local maximum *)
       Hess = D[a, {indepVars, 2}];
       localmax = Select[res,
         NegativeDefiniteMatrixQ[
           TimeConstrained[Simplify[Hess /. #], 60, Hess /. #]] &];

       If[Length[localmax] == 0,
        Throw[<|"n" -> nCopies, "w" -> w,
          "consistency" -> "a has no local maximum"|>]
       ];

       bestSol = Last[SortBy[localmax, (a /. #) &]];
       Scan[(q[#[[1]] /. q[sym_] :> sym] = #[[2]]) &, bestSol];
      ];
     ];

     (* ===== 5. Consistency Checks ===== *)
     If[!Element[N[SetPrecision[a, 50]], Reals] ||
        !Element[N[SetPrecision[c, 50]], Reals],
      Throw[<|"n" -> nCopies, "w" -> w,
        "consistency" -> "imaginary central charges"|>]
     ];

     rcharge = AssociationThread[
       ToString /@ flatMat ->
        Expand[Simplify[q[#]] & /@ flatMat]
     ];

     consistency =
      If[Min[Values[N[SetPrecision[rcharge, 50]]]] <= 0 ||
         Min[N[SetPrecision[{a, c}, 50]]] <= 0 ||
         MemberQ[Values[N[SetPrecision[rcharge, 50]]], 2.`50],
       <|"consistency" -> "non-positive R-charges"|>,
       <|"consistency" -> "consistent"|>
      ];

     If[N[SetPrecision[c, 50]] != 0 &&
        (1/2 > N[SetPrecision[a/c, 50]] ||
         N[SetPrecision[a/c, 50]] > 3/2),
      consistency = <|"consistency" -> "Hofman-Maldacena bound violated"|>
     ];

     (* ===== 6. Output ===== *)
     central = <|"a" -> a, "c" -> c|>;

     rational = If[
       (SameQ[Rational, Head[Rationalize[a]]] || IntegerQ[Rationalize[a]]) &&
       (SameQ[Rational, Head[Rationalize[c]]] || IntegerQ[Rationalize[c]]),
       <|"rational" -> Join[
         <|"a" -> ToString[InputForm[Rationalize[a]]],
           "c" -> ToString[InputForm[Rationalize[c]]]|>,
         AssociationThread[Keys[rcharge] ->
           (ToString[InputForm[Rationalize[#]]] & /@ Values[rcharge])]
       ]|>,
       <|"rational" -> ""|>
     ];

     result = Quiet[Join[
       <|"n" -> nCopies, "w" -> w|>,
       consistency, method,
       N[SetPrecision[central, 50], 30],
       N[SetPrecision[rcharge, 50], 30],
       rational
     ]];

     result

    ] (* end Catch *)
   ]; (* end AbsoluteTiming *)

   Return[Join[finalResult, <|"computation_time" -> time|>]];
  ] (* end Module *)
 ] (* end Block *)


(* Test with SU(2) + 4 fundamentals (known: a-max gives r = 1/2 for free theory) *)
(* Uncomment to test:
testResult = FindCharges[
  3,    (* dimG = 3 for SU(2) *)
  2,    (* Tadj = h^v = 2 *)
  {{"q", 1/2, 2, 4}},  (* 4 fundamentals of SU(2), T=1/2, dim=2 *)
  {}    (* no superpotential *)
];
Print[testResult];
*)
