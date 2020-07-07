::
::  mma/interface/MathLink.tm
::
::  Copyright (C) 2020 Mario Prausa
::
::  This program is free software: you can redistribute it and/or modify
::  it under the terms of the GNU General Public License as published by
::  the Free Software Foundation, either version 3 of the License, or
::  (at your option) any later version.
::
::  This program is distributed in the hope that it will be useful,
::  but WITHOUT ANY WARRANTY; without even the implied warranty of
::  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
::  GNU General Public License for more details.
::
::  You should have received a copy of the GNU General Public License
::  along with this program.  If not, see <http://www.gnu.org/licenses/>.
::

:Evaluate:	BeginPackage["IInt`"]

:Evaluate:      IIntInit::usage = "IIntInit['default_precision'] initialize IInt"
:Evaluate:      IIntCreate::usage = "IIntCreate['kernels','x0'] create IInt object"
:Evaluate:      IIntMatchEuclidean::usage = "IIntMatchEuclidean['q',a'] match coefficients in the euclidean region"
:Evaluate:      IIntMatchPhysical::usage = "IIntMatchPhysical['q',a'] match coefficients in the physical region"
:Evaluate:      IIntSave::usage = "IIntSave['filename'] save all matching coefficients"
:Evaluate:      IIntLoad::usage = "IIntLoad['filename'] load matching coefficients"
:Evaluate:      IIntEvaluate::usage = "IIntEvaluate['identifier','x'] evaluate IInt object at 'x'"
:Evaluate:      IIntSeries::usage = "IIntSeries['identifier','var','x0','order'] expand IInt object at 'x0' upto 'order'"

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
:Function:      IIntSave
:Pattern:       IIntSave[fn_String]
:Arguments:     {fn}
:ArgumentTypes: {String}
:ReturnType:    Null
:End:

:Begin:
:Function:      IIntLoad
:Pattern:       IIntLoad[fn_String]
:Arguments:     {fn}
:ArgumentTypes: {String}
:ReturnType:    Null
:End:

:Begin:
:Function:      IIntEvaluate
:Pattern:       IIntEvaluate[identifier_String,x_String]
:Arguments:     {identifier,x}
:ArgumentTypes: {String,String}
:ReturnType:    Manual
:End:

:Begin:
:Function:      IIntSeries
:Pattern:       IIntSeries[identifier_String,var_String,x_String,order_Integer]
:Arguments:     {identifier,var,x,order}
:ArgumentTypes: {String,String,String,Integer}
:ReturnType:    Manual
:End:

:Evaluate:	End[]
:Evaluate:	EndPackage[]
