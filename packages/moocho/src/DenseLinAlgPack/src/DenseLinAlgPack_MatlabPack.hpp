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

#ifndef MATLAB_PACK_H
#define MATLAB_PACK_H

#include "DenseLinAlgPack_Types.hpp"

namespace DenseLinAlgPack {
namespace MatlabPack {

/* * @name MatlabPack.
  *
  * This package contains functions that allow integration with 
  * Matlab.
  */
// @{

/* * @name Output matrices and vectors in Matlab readable format.
  *
  * These function output vectors and matrices with enought digits
  * to reproduce the exact same floating point numbers.
  */
// @{

/** \brief . */
/* * Output a DVectorSlice.
  *
  * The DVectorSlice is output in the following format:
  *
  \verbatim
  name = [ vs(1); vs(2); ... vs(n); ];
  \endverbatim
  *
  * Above #n = vs.size()# and #'# is appended to the end if #trans != BLAS_Cpp::no_trans#.
  * Also, a newline character #\n# is appended to the end after #']'#.
  */
std::ostream& out( std::ostream& o, const char* name, const DVectorSlice& vs
  , BLAS_Cpp::Transp trans = BLAS_Cpp::no_trans );

/** \brief . */
/* * Output a DMatrixSlice.
  *
  * The DMatrixSlice is output in the following format:
  *
  \verbatim
  name = [
  gms(1,1), gms(1,2), ... gms(1,n);
  gms(2,1), gms(2,2), ... gms(2,n);
  ...       ...       ... ...
  gms(m,1), gms(m,2), ... gms(m,n);
   ];
  \endverbatim
  *
  * Above #m = gms.rows()#, #n = gms.cols()# and #'# is appended to the end
  * if #trans != BLAS_Cpp::no_trans#.
  * Also, a newline character #\n# is appended to the end
  */
std::ostream& out( std::ostream& o, const char* name, const DMatrixSlice& gms
  , BLAS_Cpp::Transp trans = BLAS_Cpp::no_trans );

// @}

// @}

}	// end namespace MatlabPack
}	// end namespace DenseLinAlgPack

#endif	// MATLAB_PACK_H
