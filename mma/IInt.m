BeginPackage["IInt`"];

Install[DirectoryName[$InputFileName]<>"/interface/IInt"];

digits = 31;

arbString[x_?NumericQ] := StringReplace[ToString[N[x,digits],InputForm],RegularExpression["``?[0-9.]+"]->""];

Options[IIntInit] = {Verbose->False};
IIntInit[prec_Integer,OptionsPattern[]] := (digits = Ceiling[prec*Log[2]/Log[10]]; IIntInit[prec,Boole[OptionValue[Verbose]]];)

IIntCreate[kernels_List,{x0_?NumericQ,x_}] := IIntObj[IIntCreate[
    StringReplace[ToString[kernels,InputForm],{"{"->"","}"->""," "->"","["->"(","]"->")"}],
    arbString[x0]],x];

IIntCreate[kernels_List,x_] := IIntCreate[kernels,{1,x}];

IIntMatchEuclidean[q_?NumericQ,a_?NumericQ] := IIntMatchEuclidean[ToString[N[q]],ToString[N[a]]];
IIntMatchEuclidean[q_?NumericQ] := IIntMatchEuclidean[ToString[N[q]],"default"];

IIntMatchPhysical[q_?NumericQ,a_?NumericQ] := IIntMatchPhysical[ToString[N[q]],arbString[a]];
IIntMatchPhysical[q_?NumericQ] := IIntMatchPhysical[ToString[N[q]],"default"];

IIntMatch[{x1_?NumericQ,x2_?NumericQ,x3_?NumericQ},q_?NumericQ] := IIntMatch[arbString[x1],arbString[x2],arbString[x3],ToString[N[q]]];

IIntEvaluate[id_String,x_?NumericQ] := ToExpression[IIntEvaluate[id,arbString[x]]];

EndPackage[];
