// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
//

#ifndef LINALGPACK_IO_NAME_LOOKUPS_H
#define LINALGPACK_IO_NAME_LOOKUPS_H

#ifdef _WINDOWS

//  MS VC++ 5.0 is not performing function look for these templated operator
//  functions as it should so I will do it for the compiler.  These
//	inline functions are injected into the local namespace.  These should
//  be removed once a standards conforming compiler is available.

namespace {

// Define some boiler plate macros

#define OPEATOR_FUNCTION(OPERATOR,STREAM_TYPE,FORMAT_TYPE,OBJECT_TYPE)								\
  inline STREAM_TYPE & OPERATOR ( STREAM_TYPE & s													\
    , DenseLinAlgPack::LinAlgPackIO:: ## FORMAT_TYPE ## <DenseLinAlgPack:: ## OBJECT_TYPE ## >& bf)	\
  {																								\
    return DenseLinAlgPack:: ## OPERATOR ## (s,bf);													\
  }

#define INPUT_OPEATOR_FUNCTION(FORMAT_TYPE,OBJECT_TYPE)												\
  OPEATOR_FUNCTION( operator>> , std::istream , FORMAT_TYPE , OBJECT_TYPE )						\

#define OUTPUT_OPEATOR_FUNCTION(FORMAT_TYPE,OBJECT_TYPE)											\
  OPEATOR_FUNCTION( operator<< , std::ostream , FORMAT_TYPE , OBJECT_TYPE )						\


INPUT_OPEATOR_FUNCTION(		bound_format		,	DVector			)
OUTPUT_OPEATOR_FUNCTION(	bound_format		,	DVector			)
OUTPUT_OPEATOR_FUNCTION(	const_bound_format	,	DVector			)

INPUT_OPEATOR_FUNCTION(		bound_format		,	DVectorSlice		)
OUTPUT_OPEATOR_FUNCTION(	bound_format		,	DVectorSlice		)
OUTPUT_OPEATOR_FUNCTION(	const_bound_format	,	DVectorSlice		)

INPUT_OPEATOR_FUNCTION(		bound_format		,	DMatrix		)
OUTPUT_OPEATOR_FUNCTION(	bound_format		,	DMatrix		)
OUTPUT_OPEATOR_FUNCTION(	const_bound_format	,	DMatrix		)

INPUT_OPEATOR_FUNCTION(		bound_format		,	DMatrixSlice	)
OUTPUT_OPEATOR_FUNCTION(	bound_format		,	DMatrixSlice	)
OUTPUT_OPEATOR_FUNCTION(	const_bound_format	,	DMatrixSlice	)

#undef OPEATOR_FUNCTION
#undef INPUT_OPEATOR_FUNCTION
#undef OUTPUT_OPEATOR_FUNCTION

}	// end namespace

#endif

#endif // LINALGPACK_IO_NAME_LOOKUPS_H
