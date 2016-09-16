(* ::Package:: *)

singleRun=Function[{startEq,newInits,obfun},
Block[{newstate=First@NDSolve`Reinitialize[startEq,newInits],sol,values},
NDSolve`Iterate[newstate,tmax];
sol=NDSolve`ProcessSolutions[newstate][[All,2]];
values=Outer[#1[#2]&,sol,times];
Chop[obfun[values]]
]
];


singleRunShort=Function[{startEq,newInits,times},
Block[{newstate=First@NDSolve`Reinitialize[startEq,newInits],sol,values},
NDSolve`Iterate[newstate,Last[times]];
sol=NDSolve`ProcessSolutions[newstate][[All,2]];
values=Outer[#1[#2]&,sol,times]
]
];
