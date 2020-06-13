:Evaluate:	BeginPackage["IInt`"]

:Evaluate:      IIntInit::usage = "IIntInit['default_precision'] initialize IInt"
:Evaluate:      IIntCreate::usage = "IIntCreate['kernels','x0'] create IInt object"
:Evaluate:      IIntMatch::usage = "IIntMatch['x1','x2','x3','q'] match coefficients 'x1'->'x2'->'x3'"
:Evaluate:      IIntEvaluate::usage = "IIntEvaluate['identifier','x'] evaluate IInt object at 'x'"

:Evaluate:	Begin["Private`"]

:Begin:
:Function:      IIntInit
:Pattern:       IIntInit[iprec_Integer,iverbose_Integer]
:Arguments:     {iprec,iverbose}
:ArgumentTypes: {Integer,Integer}
:ReturnType:    Null
:End:

:Begin:
:Function:      IIntCreate
:Pattern:       IIntCreate[kernels_String,x0_String]
:Arguments:     {kernels,x0}
:ArgumentTypes: {String,String}
:ReturnType:    Manual
:End:

:Begin:
:Function:      IIntMatch
:Pattern:       IIntMatch[x1_String,x2_String,x3_String,q_String]
:Arguments:     {x1,x2,x3,q}
:ArgumentTypes: {String,String,String,String}
:ReturnType:    Null
:End:

:Begin:
:Function:      IIntEvaluate
:Pattern:       IIntEvaluate[identifier_String,x_String]
:Arguments:     {identifier,x}
:ArgumentTypes: {String,String}
:ReturnType:    Manual
:End:

:Evaluate:	End[]
:Evaluate:	EndPackage[]
