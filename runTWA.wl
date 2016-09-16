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


<<constBose2dSU3.wl


(*<<constJonathan.wl*)


(*<<constRepSU3.wl*)


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


(*Timing[
fullTWA=0;
squares=0;
Table[
newObs=singleRun[start,Flatten[discInitsOR],obsfun];
AddTo[fullTWA,newObs/runs];
AddTo[squares,newObs^2/runs];
,{rr,runs}];]*)


Timing[
fullTWA=0;
squares=0;
firstTime=First@splitTimes;
nextTimes=Drop[splitTimes,1];
Table[
stuff=singleRunShort[start,discInitsOR,firstTime];
resR=Transpose[stuff[[1]],{2,3,1}];
resI=Transpose[stuff[[2]],{2,3,1}];
(*allVars={singleRunShort[start,discInitsOR,firstTime]};*)
Table[
stuff=singleRunShort[start,Flatten[discInitsMid[First[trange]-dt,{Last[resR],Last[resI]}]],trange];
resR=Join[resR,Transpose[stuff[[1]],{2,3,1}]];
resI=Join[resI,Transpose[stuff[[2]],{2,3,1}]];
,{trange,nextTimes}];
allVars={Transpose[resR,{3,1,2}],Transpose[resI,{3,1,2}]};
newObs=Chop[obsfun[allVars]];
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
