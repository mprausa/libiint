(*
 *  mma/IInt.m
 *
 *  Copyright (C) 2020 Mario Prausa
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *)

BeginPackage["IInt`"];

Install[DirectoryName[$InputFileName]<>"/interface/IInt"];

digits = 31;

arbString[x_Real] := StringReplace[ToString[N[x,digits],InputForm],{RegularExpression["``?[0-9.]+"]->"","*^"->"e"}];
arbString[x_Integer] := ToString[x,InputForm];
arbString[x_Complex] := arbString[Re[x]]<>" + "<>arbString[Im[x]]<>"*I";
arbString[x_?NumericQ] /; Head[x] =!= Real && Head[x] =!= Complex && Head[x] =!= Integer := arbString[N[x,digits]];

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

IIntEvaluate[id_String,x_?NumericQ] := ToExpression[IIntEvaluate[id,arbString[x]]];

IIntSeries[id_String,{x_Symbol,x0_?NumericQ,order_Integer}] := ToExpression[IIntSeries[id,ToString[x],arbString[x0],order]];

EndPackage[];
