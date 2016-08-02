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


runs=1;


(* ::Subsubsection:: *)
(*vars*)


length=18;


sites=length;


suLocalDim=2;


clustSize=1;


suClustDim=suLocalDim^clustSize;


numClust=sites/clustSize;


ordPairs=Flatten[Table[Table[{ii,jj},{jj,ii,suClustDim}],{ii,suClustDim}],1];


bonds=Table[{n,Mod[n+1,length,1]},{n,length-1}];


(*bonds={{1,2}};*)


(* ::Subsubsection:: *)
(*Ham*)


localHam=PauliMatrix[3];


dis=20;


localPot=RandomReal[{-dis,dis},length];


(*crosHam1=PauliMatrix[1];
crosHam2=PauliMatrix[2];*)


crosHam1=PauliMatrix/@Range[3];
crosHam2=PauliMatrix/@Range[3];


initKet=Riffle[{0,1}&/@Range[numClust/2],{1,0}&/@Range[numClust/2]];


(* ::Subsubsection:: *)
(*Observable*)


siteObs=PauliMatrix[3];


basicObs=Table[cO[ss][#1[[1]],#1[[2]]],{ss,numClust}]&@@@Drop[ArrayRules[siteObs],-1];


multBy=#2&@@@Drop[ArrayRules[siteObs],-1];


observables=basicObs;
(*{cO[#][1,1]&/@Range[numClust]
(*Flatten[{Em[#1,#2]&@@@midPairs,El[#1,#2]&@@@lowPairs}]*)
};*)
obsfun=Function[{values},
{multBy.values,(multBy.values).Table[(-1)^n,{n,numClust}],((multBy.values).Table[(-1)^n,{n,numClust}])^2
(*values[[1]]+1/2*)(*Total[(values[[1]]\[Transpose])^2]/2+Total[Abs[values[[3]]\[Transpose]]^2]*)}
];
