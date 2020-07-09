(*
 * eval/match_physical.m
 *
 * Coypright (C) 2020 Mario Prausa
 *)

(* match all iterated integrals of cggh[2,2] along a physical path *)

<<IInt`

(* initialize IInt with a working precision of 200 bits *)
IIntInit[200,Verbose->True]

<<ggh-nc2.m

(* register all IInt's in cggh[2,2] *)
IIntCreate[cggh[2,2]];

(* match along physical path *)
IIntMatchPhysical[0.0078125];

(* save constants *)
IIntSave["physical.yaml"];






