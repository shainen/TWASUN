(* ::Package:: *)

(*tscale=10;*)


tmax=300;
steps=1000;
dt=tmax/(steps-1);
times=Range[0,tmax,dt];
split=100;
splitTimes=Partition[times,steps/split];


(*tminExp=-2;
tmaxExp=2;
tmax=10.^tmaxExp;
steps=1000;
tExps=Range[tminExp,tmaxExp,(tmaxExp-tminExp)/(steps-1)];
times=10.^#&/@tExps;*)


runs=1;


(* ::Subsubsection:: *)
(*vars*)


length=31;


sites=length^2;


suLocalDim=3;


clustSize=1;


suClustDim=suLocalDim^clustSize;


numClust=sites/clustSize;


containedSite[clust_,num_]:=clustSize(clust-1)+num


containedSites[clust_]:=Table[containedSite[clust,n],{n,clustSize}]


(*ordPairs=Flatten[Table[Table[{ii,jj},{jj,ii,suClustDim}],{ii,suClustDim}],1];*)


realPairs=Flatten[Table[Table[{ii,jj},{jj,ii,suClustDim}],{ii,suClustDim}],1];


imPairs=Flatten[Table[Table[{ii,jj},{jj,ii+1,suClustDim}],{ii,suClustDim}],1];


coToLiR[co_]:=Position[realPairs,co]


coToLiI[co_]:=Position[imPairs,co]


(*bonds=Table[{n,Mod[n+1,length,1]},{n,length-1}];*)


siteToClustNum[site_]:=Quotient[site,clustSize,1]+1


siteToClustPos[site_]:=Mod[site,clustSize,1]


addl[num_]:=Mod[num,length]


nfc[coord_]:=FromDigits[addl[coord],length]+1


bondsHoriz=Flatten[Table[{nfc[{xx,yy}],nfc[{xx,yy}+{0,1}]},{xx,0,length-1},{yy,0,length-2}],1];


bondsVert=Flatten[Table[{nfc[{xx,yy}],nfc[{xx,yy}+{1,0}]},{xx,0,length-2},{yy,0,length-1}],1];


bonds=Join[bondsHoriz,bondsVert];


extBonds=Complement[bonds,Flatten[Tuples[#,2]&/@Table[Table[containedSite[cc,ss]
,{ss,clustSize}],{cc,numClust}],1]];


intBonds=Complement[bonds,extBonds];


(*clustFromSite[ss_]:=Quotient[ss-1,clustSize]+1*)


(*clustBonds=Map[clustFromSite,extBonds,{2}];*)


(*bonds={{1,2}};*)


(* ::Subsubsection:: *)
(*Random Pot*)


norm1d=Table[PDF[NormalDistribution[0,0.5],x],{x,-(length-1)/2,(length-1)/2}];


kern=Outer[Times,norm1d,norm1d];


sqrands=RandomReal[{0,1},{2*length-1,2*length-1}]^2;


corrands=ListConvolve[kern,sqrands];


(* ::Subsubsection:: *)
(*Params*)


dis=8;


intU=24.4;


mass=86.9*1.66*10^-27;


alat=532*10^-9;


hopJ=24.8*6.626*10^-34;


harmConst=mass*alat^2/hopJ;


harmV=harmConst*Table[((2\[Pi]*54)^2 x^2+(2\[Pi]*60)^2 y^2)/2,{x,-(length-1)/2,(length-1)/2},{y,-(length-1)/2,(length-1)/2}];


(* ::Subsubsection:: *)
(*Matrices*)


bop=SparseArray[{i_,j_}/;j-i==1:>Sqrt[i],{suLocalDim,suLocalDim}];


numOp=bop\[ConjugateTranspose].bop;


sX=Sqrt[2](bop+bop\[ConjugateTranspose])/2;


sY=Sqrt[2](bop-bop\[ConjugateTranspose])/(2I);


intMat=intU/2*numOp.(numOp-IdentityMatrix[suLocalDim]);


(* ::Subsubsection:: *)
(*Ham*)


clustOp[op_,ss_]:=KroneckerProduct[IdentityMatrix[suLocalDim^(ss-1)],op,IdentityMatrix[suLocalDim^(clustSize-ss)]]


siteOp[op_,ss_]:=clustOp[op,siteToClustPos[ss]]


disConst=4*dis;


(*crosCoup[s1_,s2_]:=Abs[s1-s2]^-\[Alpha]+(length-Abs[s1-s2])^-\[Alpha]*)


(*crosCoup[s1_,s2_]:=Min[Abs[s1-s2],(length-Abs[s1-s2])]^-\[Alpha]*)


(*crosCoup[s1_,s2_]:=If[Abs[s1-s2]\[Equal]1,1,0]*)


crosCoup[s1_,s2_]:=-1


localPot=Flatten[disConst*corrands+harmV];


(*localPot=Array[h,length];*)


potHam=Table[Sum[localPot[[containedSite[cc,ss]]]
clustOp[numOp,ss]
,{ss,clustSize}],{cc,numClust}];


selfCoup=Table[
Sum[
If[MemberQ[intBonds,{containedSite[cc,s1],containedSite[cc,s2]}],
crosCoup[s1,s2](clustOp[sX,s1].clustOp[sX,s2]+clustOp[sY,s1].clustOp[sY,s2]),
0
]
,{s1,clustSize},{s2,clustSize}],{cc,numClust}];


localHam=potHam+selfCoup;


(*crosHam1=PauliMatrix[1];
crosHam2=PauliMatrix[2];*)


(* ::Input:: *)
(*(*crosHam={clustOp[PauliMatrix[#],clustSize]&/@Range[3],clustOp[PauliMatrix[#],1]&/@Range[3]};*)*)


crosHamFunc[site_]:={siteOp[sX,site],siteOp[sY,site]}


initKetSite=Flatten[Table[If[x^2+y^2<9^2&&y<=0,{0,1,0},{1,0,0}],{x,-(length-1)/2,(length-1)/2},{y,-(length-1)/2,(length-1)/2}],1];


(*initKet=Table[
Flatten[KroneckerProduct@@(List/@initKetSite[[containedSite[cc,1];;containedSite[cc,clustSize]]]),1]
,{cc,numClust}];*)


initKet=initKetSite;


(* ::Subsubsection:: *)
(*Observable*)


(*siteObs=PauliMatrix[3]*)


(*basicObs=Table[cR[ss][#,#],{ss,numClust}]&/@Range[suClustDim];*)


allObs={Table[cR[ss]@@@realPairs,{ss,numClust}]\[Transpose],Table[cI[ss]@@@imPairs,{ss,numClust}]\[Transpose]};


multBy=Table[Diagonal[clustOp[numOp,n]],{n,clustSize}];


(*multByCross=Table[Diagonal[clustOp[PauliMatrix[3],n[[1]]].clustOp[PauliMatrix[3],n[[2]]]],{n,Flatten[Table[Table[{ii, jj}, {jj, ii+1, clustSize}], {ii, clustSize}], 1]}];*)


matToVars[mat_, clVars_] := Total[(If[Length[coToLiR[#1]] == 1, 2 Re[#2] clVars[[1, coToLiR[#1][[1, 1]]]], 0]/(1 + KroneckerDelta[#1[[1]], #1[[2]]]) + If[Length[coToLiI[#1]] == 1, -2 Im[#2] clVars[[2, coToLiI[#1][[1, 1]]]], 0]) & @@@ Drop[ArrayRules[mat], -1]]


observables=allObs;
obsfun=Function[{values},
(*avNum=Flatten[(#.values&/@multBy)\[Transpose],1];*)
avNum=matToVars[numOp,values];
leftmright=Flatten[Table[If[y<=0,1,-1],{x,5},{y,-(length-1)/2,(length-1)/2}]];
numToAv=5;
imb=leftmright.avNum[[length*(length-numToAv)/2+1;;length*(length+numToAv)/2]];
{avNum,imb}
];
