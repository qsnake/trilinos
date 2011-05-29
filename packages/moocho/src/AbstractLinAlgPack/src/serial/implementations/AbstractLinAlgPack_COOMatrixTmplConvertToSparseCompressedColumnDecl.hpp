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

#ifndef COO_MATRIX_TMPL_CONVERT_TO_SPARSE_COMPRESSED_COLUMN_DECL_H
#define COO_MATRIX_TMPL_CONVERT_TO_SPARSE_COMPRESSED_COLUMN_DECL_H

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** @name {\bf Conversion to Fortran compatable sparse compressed column 
  * operations for COOMatrixTemplateInterface (Level 2,3 BLAS)}.
  *
  * See the ConvertToCSC class.
  */
//@{

/** \brief . */
template<class T_COOM>
size_type COOM_num_in_column(
    const T_COOM&						m
  , BLAS_Cpp::Transp					trans
  , size_type							col_offset
  , const IVector::value_type*		col_perm
  , size_type*						num_in_col	);

/** \brief . */
template<class T_COOM>
void COOM_insert_nonzeros(
    const T_COOM&						m
  , BLAS_Cpp::Transp					trans
  , value_type						alpha
  , size_type							row_offset
  , size_type							col_offset
  , const IVector::value_type*		row_perm
  , const IVector::value_type*		col_perm
  , size_type*						next_nz_in_col
  , FortranTypes::f_dbl_prec*			D_val
  , FortranTypes::f_int*				D_row_i			);

/** \brief . */
template<class T_COOM>
value_type COOM_insert_scaled_nonzeros(
    const T_COOM&						m
  , BLAS_Cpp::Transp					trans
  , value_type						scaled_max_ele
  , size_type							row_offset
  , size_type							col_offset
  , const IVector::value_type*		row_perm
  , const IVector::value_type*		col_perm
  , size_type*						next_nz_in_col
  , FortranTypes::f_dbl_prec*			D_val
  , FortranTypes::f_int*				D_row_i			);

//@}

} // end namespace AbstractLinAlgPack

#endif	// COO_MATRIX_TMPL_CONVERT_TO_SPARSE_COMPRESSED_COLUMN_DECL_H
