(* ::Package:: *)

(* ::Section:: *)
(*All*)


(* ::Subsubsection:: *)
(*setup*)


(*SetDirectory[NotebookDirectory[]]*)


SetDirectory[Directory[]<>"/TWASUN"];


<<randomSeed.wl


<<dynSUNRI.wl


<<ndsolve.wl


<<inits.wl


<<gausInits.wl


(*<<constBose2dSU3.wl*)


(*<<constJonathan.wl*)


<<constRepSU3.wl


(*<<constHeisAllToAllSpin1.wl*)


(* ::Subsection:: *)
(*run TWA*)


Dynamic[rr]


Timing[start=makeDSolveStart[localHam,crosHamFunc,observables];]


(*Timing[eachTWA={};
Table[
AppendTo[eachTWA,singleRun[start,Flatten[gausInitsOR],obsfun]];
,{rr,runs}];
fullTWA=Total[eachTWA]/runs;
varTWA=Variance[eachTWA];]*)


Timing[
fullTWA=0;
squares=0;
Table[
newObs=singleRun[start,Flatten[discInitsOR],obsfun];
AddTo[fullTWA,newObs/runs];
AddTo[squares,newObs^2/runs];
,{rr,runs}];]


(*Timing[fullTWA=0;
Table[
<<constRandHeisLR.wl;
start=makeDSolveStart[localHam,crosHamFunc,observables];
AddTo[fullTWA,singleRun[start,Flatten[discInitsOR],obsfun]/runs];
,{rr,runs}];]*)


mmu=MaxMemoryUsed[]/10.^6;


SetDirectory[ParentDirectory[]];


allData=fullTWA;
stError=Sqrt[squares-fullTWA^2]/Sqrt[runs];


Save["dataTWA.dat",{mmu,allData,stError}];
