(* ::Package:: *)

(*vO[ss_][xx_,yy_]:=Which[yy<xx,cO[ss][yy,xx][t]\[Conjugate],True,cO[ss][xx,yy][t]]*)


vO[ss_][xx_,yy_]:=Which[xx==yy,cR[ss][xx,yy][t],yy<xx,cR[ss][yy,xx][t]-I cI[ss][yy,xx][t],True,cR[ss][xx,yy][t]+I cI[ss][xx,yy][t]]


dotO=Function[{ss,a,b,sudim},
SparseArray[{{i_,a}->-1.0*I vO[ss][i,b]},{sudim,sudim}]+
SparseArray[{{b,i_}->1.0*I vO[ss][a,i]},{sudim,sudim}]
];


dotMatHam[ss_,ham_]:=Total[#2 dotO[ss,#1[[1]],#1[[2]],Length[ham]]&@@@Drop[ArrayRules[ham],-1]]


varHam[ss_,ham_]:=Total[#2 vO[ss][#1[[1]],#1[[2]]]&@@@Drop[ArrayRules[ham],-1]]


crossDotHam[a_,b_,cHam1_,cHam2_]:=dotMatHam[a,cHam1]varHam[b,cHam2]


makeDSolveStart = Function[{localHam, crosHamFunc, observables},
   eqsO = Table[
     matLocal = dotMatHam[ss, localHam[[ss]]];
     matCross=Sum[
Sum[
otherSite=(Complement[pairs, {inSites}][[1]]);
crosCoup[inSites,otherSite]Sum[crossDotHam[ss,siteToClustNum[otherSite],crosHamFunc[inSites][[n]], crosHamFunc[otherSite][[n]]],{n,Length[crosHamFunc[inSites]]}]
,{pairs,Cases[extBonds,{_,inSites}|{inSites,_}]}]
,{inSites, containedSites[ss]}];
     matEq = matLocal + matCross;
     {cR[ss][#1, #2]'[t] == Chop[Simplify[(matEq[[#1, #2]] + matEq[[#2, #1]])/2]] & @@@ realPairs , cI[ss][#1, #2]'[t] == Chop[Simplify[(matEq[[#1, #2]] - matEq[[#2, #1]])/(2 I)]] & @@@ imPairs}
     , {ss, numClust}];
   initsO = Table[{cR[ss][#1, #2][0] == 0 & @@@ realPairs , cI[ss][#1, #2][0] == 0 & @@@ imPairs}, {ss, numClust}];
                        First@NDSolve`ProcessEquations[Flatten[{eqsO, initsO}], observables, t, Method -> {"EquationSimplification" -> "Solve"}]
                        ];
