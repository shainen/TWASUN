(* ::Package:: *)

vO[ss_][xx_,yy_]:=Which[yy<xx,cO[ss][yy,xx][t]\[Conjugate],True,cO[ss][xx,yy][t]]


dotO=Function[{ss,a,b,sudim},
SparseArray[{{i_,a}->-I vO[ss][i,b]},{sudim,sudim}]+
SparseArray[{{b,i_}->I vO[ss][a,i]},{sudim,sudim}]
];


dotMatHam[ss_,ham_]:=Total[#2 dotO[ss,#1[[1]],#1[[2]],Length[ham]]&@@@ArrayRules[ham]]


varHam[ss_,ham_]:=Total[#2 vO[ss][#1[[1]],#1[[2]]]&@@@ArrayRules[ham]]


crossDotHam[a_,b_,cHam1_,cHam2_]:=Total[dotMatHam[a,#1]varHam[b,#2]&@@@({cHam1,cHam2}\[Transpose])]


makeDSolveStart = Function[{localHam, crosHam1, crosHam2, pairs, observables},
eqsO = Table[matLocal = localPot[[ss]] dotMatHam[ss, localHam]+Sum[crossDotHam[ss,s2, crosHam1, crosHam2],{s2,bonds[[#1,Mod[#2+1,2,1]]]&@@@Position[bonds,ss]}];cO[ss][#1, #2]'[t] == matLocal[[#1, #2]] & @@@ ordPairs, {ss, numClust}];
         initsO = Table[cO[ss][#1, #2][0] == 0 & @@@ ordPairs, {ss, numClust}];
         First@NDSolve`ProcessEquations[Flatten[{eqsO,initsO}], observables, t, Method -> {"EquationSimplification" -> "Solve"}]
         ];
