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

#ifndef MAT_VEC_COMPARE_H
#define MAT_VEC_COMPARE_H

#include <limits>
#if defined(_GNU_GXX)
#include <cmath>
#else
#include <math.h>
#endif

#include "DenseLinAlgPack_Types.hpp"
#include "TestingHelperPack_update_success.hpp"

namespace DenseLinAlgPack {

using TestingHelperPack::update_success;

/* * @name DVectorSlice and DMatrixSlice comparison functions.
  *
  * These functions compare the elements of two DVectorSlice or DMatrixSlice
  * objects.  If any of the corresponding elements does not obey
  * abs(ele1 - ele2) < sqrt(eps) then the functions return false, otherwise
  * they return true.  An exact test (bits) is not performed to allow for some round-off
  * error to occur and still equate.
  */
// @{

/** \brief . */
const value_type sqrt_eps
#if defined(_GNU_GXX)
  = std::sqrt(std::numeric_limits<value_type>::epsilon());
#elif defined(_CPQ_CXX)
  = ::sqrt(std::numeric_limits<value_type>::epsilon());
#else
  = ::sqrt(std::numeric_limits<value_type>::epsilon());
#endif

/** \brief . */
bool comp(const DVectorSlice& vs1, const DVectorSlice& vs2);

/** \brief . */
bool comp(const DVectorSlice& vs, value_type alpha);

/** \brief . */
bool comp(const DMatrixSlice& gms1, BLAS_Cpp::Transp trans1
  , const DMatrixSlice& gms2, BLAS_Cpp::Transp trans2);

/////
//bool comp(const DMatrixSlice& gms1, const DMatrixSlice& gms2);

inline
/** \brief . */
bool comp(const DMatrixSlice& gms1, const DMatrixSlice& gms2)
{
  return comp(gms1, BLAS_Cpp::no_trans, gms2, BLAS_Cpp::no_trans);
}

/** \brief . */
bool comp(const DMatrixSlice& gms1, value_type alpha);

/** \brief . */
bool comp(const DMatrixSliceTriEle& tri_gms1, const DMatrixSliceTriEle& tri_gms2);

/** \brief . */
bool comp(const DMatrixSliceTriEle& tri_gms1, value_type alpha);

/** \brief . */
bool comp_less(const DVectorSlice& vs, value_type alpha);

// @}

}	// end namespace DenseLinAlgPack

// ////////////////////////////////////
// Inline definitions

//inline
//bool DenseLinAlgPack::comp(const DMatrixSlice& gms1, const DMatrixSlice& gms2)
//{
//	return comp(gms1, BLAS_Cpp::no_trans, gms2, BLAS_Cpp::no_trans);
//}


#endif	// MAT_VEC_COMPARE_H
