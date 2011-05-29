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
// Utilities for pivoting VectorSlices and GenMatrixSlices.

#ifndef PIVOTVECMAT_H
#define PIVOTVECMAT_H

#include <stdexcept>

#include "DenseLinAlgPack_Types.hpp"

/* * @name {\bf DVector / Matrix Permutations}.
  *
  * These are functions for pivoting the elements of a vector and the 
  * rows and/or columns of a rectandular matrix.
  *
  * For DVector and DVectorSlice pivot funcitions the #perm# argument
  * gives the mapping from the new sequence to the old sequence.
  * The statement #i_old = perm(i_new)# gives the index in the
  * old vector.  After the permutation is performed the postcondition:
  \verbatim

    perm_vs(i) == vs(perm(i)), for i = 1,2,...,vs.size()

  \endverbatim
  * is true.
  *
  * For the DMatrix permutation functions the #row_perm# argument
  * gives the row permutations and the #col_perm# argument gives the column
  * permutations respectively.
  *
  * In each of these functions the permuted object can be the same
  * as the unpermuted object.
  */

// @{

// @Include: DenseLinAlgPack_IVector.hpp

// /////////////////////////////////////////////////////////////////////////////////////////
// Public Permutation functions

namespace DenseLinAlgPack {

/** \brief . */
/* * Initialize a permutation to the identiy permutation.
  *
  * Preconditions: <ul>
  * <li> #perm.size() > 0# (throw std::length_error)
  * </ul>
  *
  * Postconditions: <ul>
  * <li> #perm(i) = i#, for i = 1, 2, ... m #perm.size()#
  * </ul>
  */
void identity_perm(IVector* perm);

/** \brief . */
/* * Find the inverse permutation (inv_perm) given a permutation vector (perm).
  *
  * Postconditions: <ul>
  * <li> #inv_perm.size() == perm.size()#
  * </ul>
  */
void inv_perm(const IVector& perm, IVector* inv_perm);

/// Permute a DVectorSlice in place
void perm_ele(const IVector& perm, DVectorSlice* vs);

/// Perform y = P*x
void perm_ele(const DVectorSlice& x, const IVector& perm, DVectorSlice* y);

/// Perform x = P'*y
void inv_perm_ele(const DVectorSlice& y, const IVector& perm, DVectorSlice* x);

/// Permute a GenMatrixSlices rows
void perm_rows(const IVector& row_perm, DMatrixSlice* gms);

/// Permute a GenMatrixSlices columns
void perm_cols(const IVector& col_perm, DMatrixSlice* gms);

/// Permute a GenMatrixSlices rows and columns
void perm_rows_cols(const IVector& row_perm, const IVector& col_perm, DMatrixSlice* gms);

#ifdef TEUCHOS_DEBUG
extern bool PermVecMat_print;
#endif

} // end namespace DenseLinAlgPack

// @}

#endif // PIVOTVECMAT_H
