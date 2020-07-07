::
::  mma/interface/MathLink.tm
::
::  Copyright (C) 2020 Mario Prausa
::
::  This program is free software: you can redistribute it and/or modify
::  it under the terms of the GNU Lesser General Public License as published by
::  the Free Software Foundation, either version 3 of the License, or
::  (at your option) any later version.
::
::  This program is distributed in the hope that it will be useful,
::  but WITHOUT ANY WARRANTY; without even the implied warranty of
::  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
::  GNU Lesser General Public License for more details.
::
::  You should have received a copy of the GNU Lesser General Public License
::  along with this program.  If not, see <http://www.gnu.org/licenses/>.
::

:Evaluate:      BeginPackage["IInt`"]
:Evaluate:      Begin["Private`"]

:Begin:
:Function:      IIntInit
:Pattern:       InternalIIntInit[iprec_Integer,iverbose_Integer]
:Arguments:     {iprec,iverbose}
:ArgumentTypes: {Integer,Integer}
:ReturnType:    Null
:End:

:Begin:
:Function:      IIntCreate
:Pattern:       InternalIIntCreate[kernels_String,x0_String]
:Arguments:     {kernels,x0}
:ArgumentTypes: {String,String}
:ReturnType:    Manual
:End:

:Begin:
:Function:      IIntMatchEuclidean
:Pattern:       InternalIIntMatchEuclidean[q_String,a_String]
:Arguments:     {q,a}
:ArgumentTypes: {String,String}
:ReturnType:    Null
:End:

:Begin:
:Function:      IIntMatchPhysical
:Pattern:       InternalIIntMatchPhysical[q_String,a_String]
:Arguments:     {q,a}
:ArgumentTypes: {String,String}
:ReturnType:    Null
:End:

:Begin:
:Function:      IIntSave
:Pattern:       InternalIIntSave[fn_String]
:Arguments:     {fn}
:ArgumentTypes: {String}
:ReturnType:    Null
:End:

:Begin:
:Function:      IIntLoad
:Pattern:       InternalIIntLoad[fn_String]
:Arguments:     {fn}
:ArgumentTypes: {String}
:ReturnType:    Null
:End:

:Begin:
:Function:      IIntEvaluate
:Pattern:       InternalIIntEvaluate[identifier_String,x_String]
:Arguments:     {identifier,x}
:ArgumentTypes: {String,String}
:ReturnType:    Manual
:End:

:Begin:
:Function:      IIntSeries
:Pattern:       InternalIIntSeries[identifier_String,var_String,x_String,order_Integer]
:Arguments:     {identifier,var,x,order}
:ArgumentTypes: {String,String,String,Integer}
:ReturnType:    Manual
:End:

:Evaluate:      End[]
:Evaluate:      EndPackage[]
