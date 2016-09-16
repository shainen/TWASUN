(* ::Package:: *)

matR[i_,j_]:=(SparseArray[{{i,j}->1},{suClustDim,suClustDim}]+SparseArray[{{j,i}->1},{suClustDim,suClustDim}])/2


matI[i_,j_]:=(SparseArray[{{i,j}->-I},{suClustDim,suClustDim}]+SparseArray[{{j,i}->I},{suClustDim,suClustDim}])/2


randFromMat[mat_,init_]:=RandomChoice[Abs[init.#]^2&/@(Eigensystem[N[mat]][[2]])->(Eigensystem[N[mat]][[1]])]


(* ::Subsubsection:: *)
(*Complex*)


discRandO[s_,i_,j_]:=If[i==j,
randFromMat[matR[i,j],initKet[[s]]],
randFromMat[matR[i,j],initKet[[s]]]+I randFromMat[matI[i,j],initKet[[s]]]
]


discInitsO := Table[cO[ss][#1, #2][0] == discRandO[ss,#1,#2] & @@@ ordPairs, {ss, numClust}]


(* ::Subsubsection:: *)
(*Real*)


discRandOR[s_,i_,j_]:=randFromMat[matR[i,j],initKet[[s]]]


discRandOI[s_,i_,j_]:=randFromMat[matI[i,j],initKet[[s]]]


discInitsOR := Table[{cR[ss][#1, #2][0] == discRandOR[ss,#1,#2] & @@@ realPairs,cI[ss][#1, #2][0] == discRandOI[ss,#1,#2] & @@@ imPairs}, {ss, numClust}]


discInitsMid[st_,obs_] := Table[{cR[ss][#1, #2][st] == obs[[1,coToLiR[{#1,#2}][[1,1]],ss]] & @@@ realPairs,cI[ss][#1, #2][st] ==obs[[2,coToLiI[{#1,#2}][[1,1]],ss]] & @@@ imPairs}, {ss, numClust}]
