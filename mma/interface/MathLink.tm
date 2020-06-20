:Evaluate:	BeginPackage["IInt`"]

:Evaluate:      IIntInit::usage = "IIntInit['default_precision'] initialize IInt"
:Evaluate:      IIntCreate::usage = "IIntCreate['kernels','x0'] create IInt object"
:Evaluate:      IIntMatchEuclidean::usage = "IIntMatchEuclidean['q',a'] match coefficients in the euclidean region"
:Evaluate:      IIntMatchPhysical::usage = "IIntMatchPhysical['q',a'] match coefficients in the physical region"
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
:Function:      IIntMatchEuclidean
:Pattern:       IIntMatchEuclidean[q_String,a_String]
:Arguments:     {q,a}
:ArgumentTypes: {String,String}
:ReturnType:    Null
:End:

:Begin:
:Function:      IIntMatchPhysical
:Pattern:       IIntMatchPhysical[q_String,a_String]
:Arguments:     {q,a}
:ArgumentTypes: {String,String}
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
