(* ::Package:: *)

(* ::Subsection:: *)
(*import*)


(*import[rname_,runs_,ignore_,comb_]:=(
list=Partition[Complement[Range[0,runs-1],ignore],comb];
eachOne={};
eachVarImb={};
Do[
tempAll=0;
tempVar=0;
Do[
Get["/data/shainen/"<>rname<>"/r"<>ToString[kk]<>"/spinchain.dat"];
data={TWASpDiscSU4,TWASpGauSU4,TWASpDelta};
AddTo[tempAll,data];
AddTo[tempVar,data[[All,3]]-data[[All,2]]^2];
,{kk,rr}];
AppendTo[eachOne,tempAll/comb];
AppendTo[eachVarImb,tempVar/comb];
,{rr,list}];
all=Total[eachOne]/Length[list];
avImb=all[[All,2]];
varImb=Total[eachVarImb]/Length[list];
)*)


import[rname_,runs_,ignore_,comb_]:=(
list=Partition[Complement[Range[0,runs-1],ignore],comb];
eachOne={};
Do[
tempAll=0;
Do[
Get["/data/shainen/"<>rname<>"/r"<>ToString[kk]<>"/dataTWA.dat"];
data=allData;
AddTo[tempAll,data];
,{kk,rr}];
AppendTo[eachOne,tempAll/comb];
,{rr,list}];
avg=Total[eachOne]/Length[list];
)


importE[rname_,runs_,ignore_]:=(
list=Complement[Range[0,runs-1],ignore];
tempAll=0;
tempSq=0;
Do[
Get["/data/shainen/"<>rname<>"/r"<>ToString[rr]<>"/dataTWA.dat"];
AddTo[tempAll,allData];
AddTo[tempSq,stError^2];
,{rr,list}];
mean=tempAll/Length[list];
stEr=Sqrt[tempSq]/Length[list];
)


dir=StringSplit[ParentDirectory[],"/"][[5]];


import[dir,100,{},10];


importE[dir,100,{}];


Save["/data/shainen/"<>dir<>"compiled.dat",{avg,eachOne}];


Save["/data/shainen/"<>dir<>"_compiled.dat",{mean,stEr}];
