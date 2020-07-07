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

    IInt::usage = "IInt[kernels,{x,x0}] denotes an iterated integral\nIInt[kernels,x] is equivalent to IInt[kernels,{x,1}]";
    IIntObj::usage = "IIntObj[id,x] is an object returned by IIntCreate";

    IIntInit::usage = "IIntInit[prec] initializes IInt";
    IIntCreate::usage = "IIntCreate[kernels,x] creates IIntObj for IInt[kernels,x]\nIIntCreate[ex] creates IIntObj for all IInt's in ex";
    IIntMatchEuclidean::usage = "IIntMatchEuclidean[q,a] matches constant terms of all created IIntObj's along a euclidean path, a defaults to Infinity";
    IIntMatchPhysical::usage = "IIntMatchPhysical[q,a] matches constant terms of all created IIntObj's along a physical path, a defaults to Infinity";
    IIntSave::usage = "IIntSave[fn] saves all constant terms to a YAML file fn";
    IIntLoad::usage = "IIntLoad[fn] loads all constant terms from a YAML file fn";
    IIntEvaluate::usage = "IIntEvaluate[id,x] evaluates an IInt specified by id at x\nIIntEvaluate[ex] evaluates all IIntObj's in ex";
    IIntSeries::usage = "IIntSeries[id,{x,x0,order}] expands an IInt specified by id around x0 up to (x-x0)^order\nIIntSeries[ex] expands all IIntObj's in ex";

    Begin["Private`"];
        digits = 31;

        arbString[x_Real] := StringReplace[ToString[N[x,digits],InputForm],{RegularExpression["``?[0-9.]+"]->"","*^"->"e"}];
        arbString[x_Integer] := ToString[x,InputForm];
        arbString[x_Complex] := arbString[Re[x]]<>" + "<>arbString[Im[x]]<>"*I";
        arbString[x_?NumericQ] /; Head[x] =!= Real && Head[x] =!= Complex && Head[x] =!= Integer := arbString[N[x,digits]];

        Options[IIntInit] = {Verbose->False};
        IIntInit[prec_Integer,OptionsPattern[]] := (digits = Ceiling[prec*Log[2]/Log[10]]; InternalIIntInit[prec,Boole[OptionValue[Verbose]]];)

        IIntCreate[kernels_List,{x0_?NumericQ,x_}] := IIntObj[InternalIIntCreate[
            StringReplace[ToString[kernels,InputForm],{"{"->"","}"->""," "->"","["->"(","]"->")"}],
            arbString[x0]],x];

        IIntCreate[kernels_List,x_] := IIntCreate[kernels,{1,x}];
        IIntCreate[ex_] := (ex/.IInt->IIntCreate);

        IIntMatchEuclidean[q_?NumericQ,a_?NumericQ] := InternalIIntMatchEuclidean[ToString[N[q]],ToString[N[a]]];
        IIntMatchEuclidean[q_?NumericQ] := InternalIIntMatchEuclidean[ToString[N[q]],"default"];

        IIntMatchPhysical[q_?NumericQ,a_?NumericQ] := InternalIIntMatchPhysical[ToString[N[q]],arbString[a]];
        IIntMatchPhysical[q_?NumericQ] := InternalIIntMatchPhysical[ToString[N[q]],"default"];

        IIntSave[fn_String] := InternalIIntSave[fn];
        IIntLoad[fn_String] := InternalIIntLoad[fn];

        IIntEvaluate[id_String,x_?NumericQ] := ToExpression[InternalIIntEvaluate[id,arbString[x]]];
        IIntEvaluate[ex_] := (ex/.IIntObj->IIntEvaluate);

        IIntSeries[id_String,{x_Symbol,x0_?NumericQ,order_Integer}] := ToExpression[InternalIIntSeries[id,ToString[x],arbString[x0],order]];
        IIntSeries[ex_,{x_Symbol,x0_?NumericQ,order_Integer}] := (ex/.IIntObj[id_String,x]:>IIntSeries[id,{x,x0,order}]);
    End[];
EndPackage[];

