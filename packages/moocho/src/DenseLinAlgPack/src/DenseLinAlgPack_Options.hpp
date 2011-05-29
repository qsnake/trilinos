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
// Options for DenseLinAlgPack compilation
//

#ifndef LINALGPACK_OPTIONS_H
#define LINALGPACK_OPTIONS_H

#include "DenseLinAlgPack_extended_value_type.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_F77_wrappers.h"

#if !defined(LINALGPACK_NO_CHECKS)

/* * @name {\bf DenseLinAlgPack Options}.
  *
  * The header file DenseLinAlgPack_Options.hpp contains the defines for several macros that
  * determine how the library is built.  The user should comment out any
  * macros that her or she does not want to be defined.  The definition of
  * these macros cause the library code to assert the preconditions documented
  * for each of the member and non-member functions and throw the listed exceptions
  * if they are not satisfied.  Precondtions are supposed to be the 
  * responcibility of the client code so the user may only want to define
  * these macros during debugging for better program verification.
  * If the user checks all of the preconditions listed in this documentation for the calls
  * to all functions then the checks performed by the library are redundant.
  */
// @{

/** \brief . */
/* * If defined the library code checks to see if subscripts are in bounds for element access
  * an subregion indexing.  If the preconditions for the subscripting operations are
  * not satisfied then the listed exceptions will be thrown.
  */
#ifndef LINALGPACK_CHECK_RANGE
#define LINALGPACK_CHECK_RANGE 1
#endif

/** \brief . */
/* * If defined the library code checks to see if the sizes of rhs arguments in expressions are compatible.
  * The exception std::length_error will be thrown if rhs sizes are not compatible.
  */
#ifndef LINALGPACK_CHECK_RHS_SIZES
#define LINALGPACK_CHECK_RHS_SIZES 1
#endif

/** \brief . */
/* * If defined the library code checks to see if DVectorSlice and DMatrixSlice objects have valid constructions.
  * If they do not have valid constructions then an exception will be thrown.  The operation of these
  * checks may depend on the definition of the macro \Ref{LINALGPACK_CHECK_RANGE}.
  */
#ifndef LINALGPACK_CHECK_SLICE_SETUP
#define LINALGPACK_CHECK_SLICE_SETUP 1
#endif

#endif

namespace DenseLinAlgPack{

/// Typedef for the value type of elements that is used for the library.
typedef FortranTypes::f_dbl_prec		value_type;
/// Typedef for the index type of elements that are used by the library
typedef Teuchos::Ordinal index_type;
/// Typedef for the size type of elements that are used by the library
typedef	Teuchos::Ordinal size_type;

}

// @}

#endif // LINALGPACK_OPTIONS_H
