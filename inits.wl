(* ::Package:: *)

matR[i_,j_]:=SparseArray[{{i,j}->1,{j,i}->1},{suClustDim,suClustDim}]/2


matI[i_,j_]:=SparseArray[{{i,j}->-I,{j,i}->I},{suClustDim,suClustDim}]/2


discRandO[s_,i_,j_]:=RandomChoice[Abs[initKet[[s]].#]^2&/@(Eigensystem[N[matR[i,j]]][[2]])->(Eigensystem[N[matR[i,j]]][[1]])]+I RandomChoice[Abs[initKet[[s]].#]^2&/@(Eigensystem[N[matI[i,j]]][[2]])->(Eigensystem[N[matI[i,j]]][[1]])]


discInitsO := Table[cO[ss][#1, #2][0] == discRandO[ss,#1,#2] & @@@ ordPairs, {ss, numClust}]
