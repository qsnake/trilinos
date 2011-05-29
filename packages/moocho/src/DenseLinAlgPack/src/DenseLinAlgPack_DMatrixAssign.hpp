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
// These are assignment functions that used to be in GenMatrixOp but because of
// some circular dependency problems their declarations where moved into a seprate
// file.

#ifndef GEN_MATRIX_ASSIGN_H
#define GEN_MATRIX_ASSIGN_H

#include "DenseLinAlgPack_Types.hpp"

namespace DenseLinAlgPack {

// ///////////////////////////////////////////////////////////////////////////////////
/* * @name {\bf DMatrix/DMatrixSlice Assignment Fucntions}.
  *
  * These are functions for assigning one matrix to another.  The rhs matrices
  * can be transposed (op(rhs) = rhs) or non-transposed (op(rhs) = trans(rhs).
  * The assignment operators for
  * DMatrix and DMatrixSlice use these functions to implement their
  * assignment operator functions.  The assignmet functions for triangular
  * matrices are ment to be called from the assignment operator functions
  * for wrapper classes for triangular (upper and lower) and sysmetric 
  * (upper and lower storage).
  *
  * These assignment functions check for aliasing and assignment to self
  * and create any needed temporaries if needed.
  */
// @{

/// gm_lhs = alpha (elementwise)
void assign(DMatrix* gm_lhs, value_type alpha);

/// gm_lhs = op(gms_rhs)
void assign(DMatrix* gm_lhs, const DMatrixSlice& gms_rhs, BLAS_Cpp::Transp trans_rhs);

/// gms_lhs = alpha (elementwise)
void assign(DMatrixSlice* gms_lhs, value_type alpha);

/// gms_lhs = op(gms_rhs)
void assign(DMatrixSlice* gms_lhs, const DMatrixSlice& gms_rhs, BLAS_Cpp::Transp trans_rhs);

/// tri_ele_gms_lhs = alpha (elementwise)
void assign(DMatrixSliceTriEle* tri_gms_lhs, value_type alpha);

/// tri_ele_gms_lhs = tri_ele_gms_rhs
void assign(DMatrixSliceTriEle* tri_ele_gms_lhs, const DMatrixSliceTriEle& tri_ele_gms_rhs);

//		end Assignment Fucntions
// @}

}	// end namespace DenseLinAlgPack

#endif	// GEN_MATRIX_ASSIGN_H
