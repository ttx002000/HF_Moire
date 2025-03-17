(* ::Package:: *)

numKernels= $ProcessorCount;


LaunchKernels[numKernels];
myPrint[x_]:=Print[NumberForm[x,{20,1}]];
Print["Available Kernels: ",Length[Kernels[]]];
Print["scriptcommandline",$ScriptCommandLine];
params=Rest[$ScriptCommandLine];
numericParams=ToExpression/@params;
U=numericParams[[1]];
J=numericParams[[2]];
Dfield=numericParams[[3]];
perturbtype=IntegerPart[numericParams[[4]]];
filepos=IntegerPart[numericParams[[5]]];
Print["U ",U];
Print["J ",J];
Print["uD ",Dfield];
Print["filepos ",filepos];
scratchDir=Environment["SCRATCH"];
outputfilename=ToString[myPrint[numericParams[[1]]]]<>"U"<>ToString[myPrint[numericParams[[2]]]]<>"J"<>ToString[myPrint[numericParams[[3]]]]<>"uD"<>ToString[myPrint[numericParams[[4]]]]<>"perturb.csv";
inputpath=FileNameJoin[{scratchDir, "stoner_RMG", "data_output"<>ToString[filepos],"stoner_perturb"<>ToString[perturbtype]<>"_input","uD"<>ToString[IntegerPart[numericParams[[3]]]]<>".csv"}];
outputpath=FileNameJoin[{scratchDir, "stoner_RMG", "data_output"<>ToString[filepos],outputfilename}];


Print["inputpath",inputpath];
Print["outputpath",outputpath]
data=Import[inputpath];
data=Rest[data];
Print["here is some of my input data",data[[1]]];
Print["here is some of my input data",data[[2]]];
uniqueData=DeleteDuplicatesBy[data,First];
x=uniqueData[[All,1]];
y=uniqueData[[All,2]];
interp=Interpolation[uniqueData,InterpolationOrder->2];
totalE[n1_,n2_,n3_,n4_,U_,J_]:=interp[n1]+interp[n2]+interp[n3]+interp[n4]-U*(n1^2+n2^2+n3^2+n4^2)+J*(n1+n2)*(n3+n4);
totaldensityrange=Range[-0.05,0.05,0.00025]   
minValues=ConstantArray[Null,Length[totaldensityrange]];
n1Values=ConstantArray[Null,Length[totaldensityrange]];   
n2Values=ConstantArray[Null,Length[totaldensityrange]];  
n3Values=ConstantArray[Null,Length[totaldensityrange]];   
n4Values=ConstantArray[Null,Length[totaldensityrange]];
SetSharedVariable[minValues,n1Values,n2Values,n3Values,n4Values];

Do[
Print["Processing density: ",totaldensityrange[[index]]];
result=NMinimize[
{totalE[n1,n2,n3,n4,U,J],Min[x]<=n1<=Max[x],
Min[x]<=n2<=Max[x],Min[x]<=n3<=Max[x],Min[x]<=n4<=Max[x],n1+n2+n3+n4==totaldensityrange[[index]],n1<=n2<=n3<=n4},{n1,n2,n3,n4},MaxIterations->1000,Method->"RandomSearch"];
Print["result ",result];
minValues[[index]]=result[[1]];    (*Store the minimum function value*)
n1Values[[index]]=(n1/. result[[2]]);
n2Values[[index]]=(n2/. result[[2]]); 
n3Values[[index]]=(n3/. result[[2]]); 
n4Values[[index]]=(n4/. result[[2]]); 
,{index,1,Length[totaldensityrange]}
];
savedata=Transpose[{n1Values,n2Values,n3Values,n4Values,minValues}];



Export[outputpath,savedata]

