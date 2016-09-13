(* ::Package:: *)

genMultNormRands=Function[{meanVec,covMat,numRands},
Block[{covnrot=Chop[Eigensystem[N[covMat]]],rands,rotRands},
rands=If[#==0.,Table[0,{numRands}],RandomVariate[NormalDistribution[0,#],numRands]]&/@Sqrt[covnrot[[1]]];
rotRands=covnrot[[2]]\[Transpose].rands+meanVec;
rotRands
]
];


gausRands=Function[{varMats,initKet,numRands},
Block[{means=initKet.#.initKet&/@varMats,covMat=Outer[N[initKet.(#1.#2+#2.#1).initKet/2-initKet.#1.initKet initKet.#2.initKet]&,varMats,varMats,1]},
genMultNormRands[means,covMat,numRands]
]
];


matR[i_,j_]:=(SparseArray[{{i,j}->1},{suClustDim,suClustDim}]+SparseArray[{{j,i}->1},{suClustDim,suClustDim}])/2


matI[i_,j_]:=(SparseArray[{{i,j}->-I},{suClustDim,suClustDim}]+SparseArray[{{j,i}->I},{suClustDim,suClustDim}])/2


gausInitsOR := Table[
matVars=Join[matR@@@realPairs,matI@@@imPairs];
siteRands=Flatten[gausRands[matVars,initKet[[ss]],1]];
{Table[cR[ss][#1, #2][0] &@@realPairs[[pp]]==siteRands[[pp]],{pp,Length[realPairs]}],Table[cI[ss][#1, #2][0] &@@imPairs[[pp]]== siteRands[[Length[realPairs]+pp]],{pp,Length[imPairs]}]}, {ss, numClust}]


meanInitsOR := Table[
matVars=Join[matR@@@realPairs,matI@@@imPairs];
siteRands=initKet[[ss]].#.initKet[[ss]]&/@matVars;
{Table[cR[ss][#1, #2][0] &@@realPairs[[pp]]==siteRands[[pp]],{pp,Length[realPairs]}],Table[cI[ss][#1, #2][0] &@@imPairs[[pp]]== siteRands[[Length[realPairs]+pp]],{pp,Length[imPairs]}]}, {ss, numClust}]
