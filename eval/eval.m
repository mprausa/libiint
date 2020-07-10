(*
 * eval/eval.m
 *
 * Coypright (C) 2020 Mario Prausa
 *)

(* This is a demonstration how to evaluate cggh[2,2] numerically via libiint *)

If[!FileExistsQ["ggh-nc2.m"],
    Print["ggh-nc2.m missing!"];
    Print["copy or symlink ggh-nc2.m into the working directory."];
    Exit[1];
];

(* working precision *)
digits = 20;

(* points *)
xpoints = {-0.00112,0.0125,-0.75};

(* small imaginary part *)
ieps = I*1.`100*^-70;

<<IInt`
<<PolyLogTools`
Print[];

(* initialize IInt *)
IIntInit[Ceiling[digits*Log[10]/Log[2]]];

(* functions to convert Plog's to GPL's *)
PlogConvert[b_List,a_List,x_] /; a =!= {} && Union[Head/@a] === {List} := Module[{sq,a1},
    sq = Sqrt[a[[1,1]]];
    a1 = Rest[a];
    Which[
        a[[1,2]] === 1,
            PlogConvert[Append[b,sq],a1,x] + PlogConvert[Append[b,-sq],a1,x],
        a[[1,2]] === 0,
            PlogConvert[Append[b,sq],a1,x]/sq - PlogConvert[Append[b,-sq],a1,x]/sq,
        True,
            $Failed
    ]
];

Plog2GPL[ex_] := (ex/.Plog[a_List,x_]:>PlogConvert[{},a,x]/.PlogConvert[a_List,{},x_]->GPL[a,x]);

(* define psi functions *)
arg1 = t*(4+t)^5/(t^2+8*t+20)/(t^2+6*t+4)^2
psi1 = 2^(1/2)*Pi*Hypergeometric2F1[1/4,3/4,1,arg1];
psi2 = I*Pi*Hypergeometric2F1[1/4,3/4,1,1-arg1];

{dpsi1,dpsi2} = Simplify[D[{psi1,psi2},t]];
{d2psi1,d2psi2} = Simplify[D[{dpsi1,dpsi2},t]];

{psi1,psi2,dpsi1,dpsi2,d2psi1,d2psi2,arg1} =
    {psi1,psi2,dpsi1,dpsi2,d2psi1,d2psi2,arg1}/.t->(1 - 9*x - Sqrt[1 - 18*x + x^2])/(2*x);

psi[x0_?NumericQ] := Which[
        Re[x0] > -9+4*Sqrt[5] && Abs[x0] < .1 && (Re[arg1/.x->x0] <= 1 || Im[arg1/.x->x0] > 0),
            psi1/.x->x0,
        True,
            (psi1 + 2*psi2)/.x->x0
      ];

dpsi[x0_?NumericQ] := Which[
        Re[x0] > -9+4*Sqrt[5] && Abs[x0] < .1 && (Re[arg1/.x->x0] <= 1 || Im[arg1/.x->x0] > 0),
            dpsi1/.x->x0,
        True,
            (dpsi1 + 2*dpsi2)/.x->x0
      ];

d2psi[x0_?NumericQ] := Which[
        Re[x0] > -9+4*Sqrt[5] && Abs[x0] < .1 && (Re[arg1/.x->x0] <= 1 || Im[arg1/.x->x0] > 0),
            d2psi1/.x->x0,
        True,
            (d2psi1 + 2*d2psi2)/.x->x0
      ];


(* load result *)
<<ggh-nc2.m

(* create IIntObj's *)
ex = IIntCreate[cggh[2,2]];

(* some replacements *)
ex = ex/.{
    r1->E^(I*Pi/3),
    r2->E^(-I*Pi/3),
    phi->GoldenRatio,
    Lmu-> I*Pi + 2*GPL[{1},x] - GPL[{0},x],  (* mu^2 = mH^2 *)
    psi->psi[x],
    dpsi->dpsi[x],
    d2psi->d2psi[x]
};

(* convert Plog's to GPL's *)
ex = Plog2GPL[ex];

(* convert to PolyLogTools notation *)
ex = ex/.GPL[{a__},x_]->G[a,x]

(* load constants *)
IIntLoad["euclidean.yaml"];
IIntLoad["physical.yaml"];

eval[x1_] := Module[{v,t1,t2},
    v = ex/.x->(x1+ieps);

    {t1,v} = AbsoluteTiming[IIntEvaluate[v]];                               (* evaluate IInt's *)
    {t2,v} = AbsoluteTiming[Ginsh[N[v,digits],{},PrecisionGoal->digits]];   (* evaluate GPL's *)

    Print[x1," => ",v," (IInt: ",t1," sec., GPL: ",t2," sec.)"];
    Return[v];
];

eval/@xpoints;
