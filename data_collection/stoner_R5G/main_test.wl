(* ::Package:: *)

numKernels=ToExpression[Environment["NUM_KERNELS"],InputForm,1];

(*Launch that many kernels*)
LaunchKernels[numKernels];
(*Print available kernels*)
Print["Available Kernels: ",Length[Kernels[]]];

(*Example parallel computation*)
result=ParallelTable[i^2,{i,1,1000}];
params=Rest[$ScriptCommandLine];
numericParams=ToExpression/@params;
U=numericParams[[1]];
J=numericParams[[2]];
Dfield=numericParams[[3]]
filepos=IntegerPart[numericParams[[4]]]
Print["U ",U];
Print["J ",J];
Print["filepos ",filepos];
scratchDir=Environment["SCRATCH"]
inputpath=scratchDir<>"stonerRMG\\data_output"<>ToString[filepos]<>"\\random_input\\input.csv"
outputpath=scratchDir<>"stonerRMG\\data_output"<>ToString[filepos]<>"\\random_output\\output.csv"
data=Import[inputpath]
data=Rest[data];
Print["inputdata",data]
uniqueData=DeleteDuplicatesBy[data,First];
dataoutput1=Transpose[{{U,J},{Dfield,filepos}}];

Export[outputpath,dataoutput1]


