(* ::Package:: *)

(*tscale=10;*)


(*tmax=4 tscale;
steps=500;
times=Range[0,tmax,tmax/(steps-1)];*)


tminExp=-2;
tmaxExp=2;
tmax=10.^tmaxExp;
steps=1000;
tExps=Range[tminExp,tmaxExp,(tmaxExp-tminExp)/(steps-1)];
times=10.^#&/@tExps;


runs=100;


(* ::Subsubsection:: *)
(*vars*)


length=18;


sites=length;


suLocalDim=2;


clustSize=2;


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


bonds=Flatten[Table[Table[{n,m},{m,n+1,length}],{n,length}],1];


(*bonds=Complement[Table[{n,Mod[n+1,length,1]},{n,length}],{{1,2}}];*)


extBonds=Complement[bonds,Flatten[Tuples[#,2]&/@Table[Table[containedSite[cc,ss]
,{ss,clustSize}],{cc,numClust}],1]];


intBonds=Complement[bonds,extBonds];


(*clustFromSite[ss_]:=Quotient[ss-1,clustSize]+1*)


(*clustBonds=Map[clustFromSite,extBonds,{2}];*)


(*bonds={{1,2}};*)


(* ::Subsubsection:: *)
(*Ham*)


clustOp[op_,ss_]:=KroneckerProduct[IdentityMatrix[suLocalDim^(ss-1)],op,IdentityMatrix[suLocalDim^(clustSize-ss)]]


siteOp[op_,ss_]:=clustOp[op,siteToClustPos[ss]]


dis=20;


\[Alpha]=2.5;


(*crosCoup[s1_,s2_]:=Abs[s1-s2]^-\[Alpha]+(length-Abs[s1-s2])^-\[Alpha]*)


crosCoup[s1_,s2_]:=Min[Abs[s1-s2],(length-Abs[s1-s2])]^-\[Alpha]


(*crosCoup[s1_,s2_]:=If[Abs[s1-s2]\[Equal]1,1,0]*)


localPot=RandomReal[{-dis,dis},length];


potHam=Table[Sum[localPot[[containedSite[cc,ss]]]
clustOp[PauliMatrix[3],ss]
,{ss,clustSize}],{cc,numClust}];


selfCoup=Table[
Sum[
If[MemberQ[intBonds,{containedSite[cc,s1],containedSite[cc,s2]}],
crosCoup[s1,s2]Sum[clustOp[PauliMatrix[ii],s1].clustOp[PauliMatrix[ii],s2],{ii,3}],
0
]
,{s1,clustSize},{s2,clustSize}],{cc,numClust}];


localHam=potHam+selfCoup;


(*crosHam1=PauliMatrix[1];
crosHam2=PauliMatrix[2];*)


(* ::Input:: *)
(*(*crosHam={clustOp[PauliMatrix[#],clustSize]&/@Range[3],clustOp[PauliMatrix[#],1]&/@Range[3]};*)*)


crosHamFunc[site_]:=siteOp[PauliMatrix[#],site]&/@Range[3]


initKetSite=Riffle[{0,1}&/@Range[sites/2],{1,0}&/@Range[sites/2]];


initKet=Table[
Flatten[KroneckerProduct@@(List/@initKetSite[[containedSite[cc,1];;containedSite[cc,clustSize]]]),1]
,{cc,numClust}];


(* ::Subsubsection:: *)
(*Observable*)


(*siteObs=PauliMatrix[3]*)


basicObs=Table[cR[ss][#,#],{ss,numClust}]&/@Range[suClustDim];


multBy=Table[Diagonal[clustOp[PauliMatrix[3],n]],{n,clustSize}];


multByCross=Table[Diagonal[clustOp[PauliMatrix[3],n[[1]]].clustOp[PauliMatrix[3],n[[2]]]],{n,Flatten[Table[Table[{ii, jj}, {jj, ii+1, clustSize}], {ii, clustSize}], 1]}];


observables=basicObs;
obsfun=Function[{values},
sigZ=Flatten[(#.values&/@multBy)\[Transpose],1];
sigCross=#.values&/@multByCross;
spinSum=Table[(-1)^n,{n,sites}].sigZ;
spinSumMid=Sum[(-1)^n sigZ[[n]],{n,2,sites,3}];
{sigZ,spinSum,spinSum^2,sigCross,spinSumMid,spinSumMid^2}
];
