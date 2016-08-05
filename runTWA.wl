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


<<constRandHeisR.wl


(* ::Subsection:: *)
(*run TWA*)


Dynamic[rr]


start=makeDSolveStart[localHam,crosHam,observables];


Timing[eachTWA={};
Table[
AppendTo[eachTWA,singleRun[start,Flatten[discInitsOR],obsfun]];
,{rr,runs}];
fullTWA=Total[eachTWA]/runs;]


mmu=MaxMemoryUsed[]/10.^6;


SetDirectory[ParentDirectory[]];


allData=fullTWA;


Save["dataTWA.dat",{mmu,allData}];
