(* ::Package:: *)

tscale=10;


tmax=4 tscale;
steps=500;
times=Range[0,tmax,tmax/(steps-1)];


runs=1;


(* ::Subsubsection:: *)
(*vars*)


length=2;


sites=length;


suLocalDim=2;


clustSize=1;


suClustDim=suLocalDim^clustSize;


numClust=sites/clustSize;


ordPairs=Flatten[Table[Table[{ii,jj},{jj,ii,suClustDim}],{ii,suClustDim}],1];


bonds=Table[{n,Mod[n+1,length,1]},{n,length}];


bonds={{1,2}};


(* ::Subsubsection:: *)
(*Ham*)


localHam=PauliMatrix[3];


(*crosHam1=PauliMatrix[1];
crosHam2=PauliMatrix[2];*)


crosHam1=PauliMatrix/@Range[3];
crosHam2=PauliMatrix/@Range[3];


initKet={1,0}&/@Range[numClust];
