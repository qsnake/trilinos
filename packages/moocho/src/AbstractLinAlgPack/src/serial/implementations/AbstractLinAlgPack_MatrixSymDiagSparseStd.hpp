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

#ifndef SPARSE_LINALG_PACK_MATRIX_DIAGONAL_SPARSE_STD_H
#define SPARSE_LINALG_PACK_MATRIX_DIAGONAL_SPARSE_STD_H

#include "AbstractLinAlgPack_MatrixSymDiagSparse.hpp"
#include "AbstractLinAlgPack_SpVectorClass.hpp"

namespace AbstractLinAlgPack {

/** \brief Concrete subclass for a serial symmetric diagonal matrix with many zeros on the diagonal.
  *
  * The underlying diagonal vector is sorted and determines the dimensions of the
  * matrix.
  *
  * The default constructor, copy constructor are allowed.
  */
class MatrixSymDiagSparseStd: virtual public MatrixSymDiagSparse {
public:

  /** @name Constructors/initializes */
  //@{

  /// Construct uninitialized
  MatrixSymDiagSparseStd()
  {}

  /** \brief Construct the diagonal.
    */
  MatrixSymDiagSparseStd( const SpVectorSlice& diag );

  /** \brief Reinitialize the diagonal.
    */
  void initialize( const SpVectorSlice& diag );

  //@}

  /** @name Overridden from MatrixOp */
  //@{

  /** \brief . */
  MatrixOp& operator=(const MatrixOp& m);

  //@}

  /** Overridden from MatrixDiagonalSparse */
  //@{

  /** \brief . */
  const SpVectorSlice diag() const;

  //@}

private:
  
  SpVector	diag_;

};	// end class MatrixDiagonalSparse

}	// end namespace AbstractLinAlgPack

#endif	// SPARSE_LINALG_PACK_MATRIX_DIAGONAL_SPARSE_STD_H
