(* ::Package:: *)

(*tscale=10;*)


tmax=100;
steps=1000;
times=Range[0,tmax,tmax/(steps-1)];


(*tminExp=-2;
tmaxExp=2;
tmax=10.^tmaxExp;
steps=1000;
tExps=Range[tminExp,tmaxExp,(tmaxExp-tminExp)/(steps-1)];
times=10.^#&/@tExps;*)


runs=100;


(* ::Subsubsection:: *)
(*vars*)


length=2;


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


(*bonds=Table[{n,Mod[n+1,length,1]},{n,length-1}];*)


siteToClustNum[site_]:=Quotient[site,clustSize,1]+1


siteToClustPos[site_]:=Mod[site,clustSize,1]


addl[num_]:=Mod[num,length]


nfc[coord_]:=FromDigits[addl[coord],length]+1


bondsHoriz=Flatten[Table[{nfc[{xx,yy}],nfc[{xx,yy}+{0,1}]},{xx,0,length-1},{yy,0,length-1}],1];


bondsVert=Flatten[Table[{nfc[{xx,yy}],nfc[{xx,yy}+{1,0}]},{xx,0,length-1},{yy,0,length-1}],1];


bonds=Join[bondsHoriz,bondsVert];


extBonds=Complement[bonds,Flatten[Tuples[#,2]&/@Table[Table[containedSite[cc,ss]
,{ss,clustSize}],{cc,numClust}],1]];


intBonds=Complement[bonds,extBonds];


(*clustFromSite[ss_]:=Quotient[ss-1,clustSize]+1*)


(*clustBonds=Map[clustFromSite,extBonds,{2}];*)


(*bonds={{1,2}};*)


(* ::Subsubsection:: *)
(*Params*)


intU=1;


hopJ=-0.1;


(* ::Subsubsection:: *)
(*Matrices*)


matSx=SparseArray[{{0,1/Sqrt[2],0},{1/Sqrt[2],0,1/Sqrt[2]},{0,1/Sqrt[2],0}}];


matSy=SparseArray[I{{0,-1/Sqrt[2],0},{1/Sqrt[2],0,-1/Sqrt[2]},{0,1/Sqrt[2],0}}];


matSz=SparseArray[{{1,0,0},{0,0,0},{0,0,-1}}];


intMat=intU/2*(matSz.matSz-IdentityMatrix[3]);


hopMats=Sqrt[2]{matSx,matSy};


(* ::Subsubsection:: *)
(*Ham*)


clustOp[op_,ss_]:=KroneckerProduct[IdentityMatrix[suLocalDim^(ss-1)],op,IdentityMatrix[suLocalDim^(clustSize-ss)]]


siteOp[op_,ss_]:=clustOp[op,siteToClustPos[ss]]


(*disConst=4*dis;*)


(*crosCoup[s1_,s2_]:=Abs[s1-s2]^-\[Alpha]+(length-Abs[s1-s2])^-\[Alpha]*)


(*crosCoup[s1_,s2_]:=Min[Abs[s1-s2],(length-Abs[s1-s2])]^-\[Alpha]*)


(*crosCoup[s1_,s2_]:=If[Abs[s1-s2]\[Equal]1,1,0]*)


crosCoup[s1_,s2_]:=hopJ;


(*localPot=Array[h,length];*)


potHam=Table[Sum[
clustOp[intMat,ss]
,{ss,clustSize}],{cc,numClust}];


selfCoup=Table[
Sum[
If[MemberQ[intBonds,{containedSite[cc,s1],containedSite[cc,s2]}],
crosCoup[s1,s2]Sum[clustOp[mats,s1].clustOp[mats,s2],{mats,hopMats}],
0
]
,{s1,clustSize},{s2,clustSize}],{cc,numClust}];


localHam=potHam+selfCoup;


(*crosHam1=PauliMatrix[1];
crosHam2=PauliMatrix[2];*)


(* ::Input:: *)
(*(*crosHam={clustOp[PauliMatrix[#],clustSize]&/@Range[3],clustOp[PauliMatrix[#],1]&/@Range[3]};*)*)


crosHamFunc[site_]:=Table[siteOp[mats,site],{mats,hopMats}]


initKetSite=Table[{0,1,0},{sites}];


(*initKet=Table[
Flatten[KroneckerProduct@@(List/@initKetSite[[containedSite[cc,1];;containedSite[cc,clustSize]]]),1]
,{cc,numClust}];*)


initKet=initKetSite;


(* ::Subsubsection:: *)
(*Observable*)


(*siteObs=PauliMatrix[3]*)


basicObs=Table[{cR[ss][1,2],cR[ss][2,3],cR[ss][1,2],cR[ss][2,3],cR[ss][1,1],cR[ss][3,3]},{ss,numClust}]\[Transpose];


(*multBy=Table[Diagonal[clustOp[numOp,n]],{n,clustSize}];*)


(*multByCross=Table[Diagonal[clustOp[PauliMatrix[3],n[[1]]].clustOp[PauliMatrix[3],n[[2]]]],{n,Flatten[Table[Table[{ii, jj}, {jj, ii+1, clustSize}], {ii, clustSize}], 1]}];*)


observables=basicObs;
obsfun=Function[{values},
avSx=Sqrt[2](values[[1]]+values[[2]]);
avSy=Sqrt[2](values[[3]]+values[[4]]);
avSz=(values[[5]]-values[[6]]);
sumSq=((Total[avSx])^2+(Total[avSy])^2-(Total[avSx^2]+Total[avSy^2]))/sites^2;
{avSx,avSy,avSz,sumSq}
];
