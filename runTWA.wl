(* ::Package:: *)

(* ::Section:: *)
(*All*)


(* ::Subsubsection:: *)
(*setup*)


(*SetDirectory[NotebookDirectory[]]*)


SetDirectory[Directory[]<>"/TWASUN"];


<<randomSeed.wl


<<dynSUNR.wl


<<ndsolve.wl


<<inits.wl


<<gausInits.wl


<<constRandHeisLR.wl


(*<<constHeisAllToAllSpin1.wl*)


(* ::Subsection:: *)
(*run TWA*)


Dynamic[rr]


Timing[start=makeDSolveStart[localHam,crosHamFunc,observables];]


(*Timing[eachTWA={};
Table[
AppendTo[eachTWA,singleRun[start,Flatten[discInitsOR],obsfun]];
,{rr,runs}];
fullTWA=Total[eachTWA]/runs;]*)


Timing[fullTWA=0;
Table[
AddTo[fullTWA,singleRun[start,Flatten[discInitsOR],obsfun]/runs];
,{rr,runs}];]


mmu=MaxMemoryUsed[]/10.^6;


SetDirectory[ParentDirectory[]];


allData=fullTWA;


Save["dataTWA.dat",{mmu,allData}];
