(* ::Package:: *)

Print["I am doing something"]
fullCommandLine = $CommandLine;

(* \:6253\:5370\:5b8c\:6574\:7684\:547d\:4ee4\:884c *)
Print["Full command line: ", fullCommandLine];


If[Length[fullCommandLine] >= 4,
    args = Drop[fullCommandLine, 3]; (* \:53bb\:6389\:524d 3 \:4e2a\:5143\:7d20\:ff08Mathematica \:8def\:5f84\:3001-script\:3001\:811a\:672c\:6587\:4ef6\:540d\:ff09 *)
    Print["Arguments: ", args],
    Print["No arguments passed."]
]
