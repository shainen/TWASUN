(* ::Package:: *)

(* ::Section:: *)
(*All*)


(* ::Subsubsection:: *)
(*setup*)


(*SetDirectory[NotebookDirectory[]]*)


SetDirectory[Directory[]<>"/TWASUN"];


<<randomSeed.wl


<<dynSUN.wl


<<ndsolve.wl


<<inits.wl


<<constRandHeis.wl


(* ::Subsection:: *)
(*run TWA*)


Dynamic[rr]


(*observables={cO[#][1,1]&/@Range[numClust]
(*Flatten[{Em[#1,#2]&@@@midPairs,El[#1,#2]&@@@lowPairs}]*)
};
obsfun=Function[{values},
{values[[1]]+1/2(*Total[(values[[1]]\[Transpose])^2]/2+Total[Abs[values[[3]]\[Transpose]]^2]*)}
];*)


(*start=makeDSolveStart[observables];*)


start=makeDSolveStart[localHam,crosHam1,crosHam2,bonds,observables];


eachTWA={};
Table[
AppendTo[eachTWA,singleRun[start,Flatten[discInitsO],obsfun]];
,{rr,runs}];
fullTWA=Total[eachTWA]/runs;


mmu=MaxMemoryUsed[]/10.^6;


SetDirectory[ParentDirectory[]];


allData=fullTWA;


Save["dataTWA.dat",{mmu,allData}];
