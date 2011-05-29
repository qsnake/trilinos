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

#ifndef INV_CHOL_UPDATE_H
#define INV_CHOL_UPDATE_H

#include "DenseLinAlgPack_Types.hpp"

namespace DenseLinAlgPack {

/** \brief . */
/* * Perform an update of the Cholesky factor of a symetric matrix (SM = L * L^T)
  * and while maintaining the triangular structure.
  *
  * This function permforms an update of L^T or L^-1 depending on whether
  * we are updating the cholesky factor or its inverse.  This function
  * implements algorithm A3.4.1 in Dennis and Schabel.  The update is:
  * (J_new^T = L_old^T + u * v^T) where J_new is rotated back to triangular form.
  *
  * Preconditions: <ul>
  * <li> UpTriM.rows() == UpTriM.cols() == u.size() == v.size() (throw std::length_error)
  * </ul>
  */
void update_chol_factor(DMatrixSlice* UpTriM, DVectorSlice* u
  , const DVectorSlice& v);

/** \brief . */
/* * Perform a jacobi rotation or a matrix about row i.
  *
  * Preconditions: <ul>
  * <li> UpTriM.rows() == UpTriM.cols() (throw std::length_error)
  * </ul>
  */
void jacobi_rotate(DMatrixSlice* UpTriM, size_type row_i, value_type alpha
  , value_type beta); 

}	// end namespace DenseLinAlgPack

#endif	// INV_CHOL_UPDATE_H
