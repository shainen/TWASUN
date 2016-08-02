(* ::Package:: *)

singleRun=Function[{startEq,newInits,obfun},
Block[{newstate=First@NDSolve`Reinitialize[startEq,newInits],sol,values},
NDSolve`Iterate[newstate,tmax];
sol=NDSolve`ProcessSolutions[newstate][[All,2]];
values=(Outer[Through[#1[#2]]&,sol,times,1]);
Chop[obfun[values]]
]
];
