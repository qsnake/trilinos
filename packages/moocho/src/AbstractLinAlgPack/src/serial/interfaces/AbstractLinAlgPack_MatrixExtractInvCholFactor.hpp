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

#ifndef MATRIX_EXTRACT_INV_CHOL_FACTOR_H
#define MATRIX_EXTRACT_INV_CHOL_FACTOR_H

#include "AbstractLinAlgPack_Types.hpp"
#include "AbstractLinAlgPack_MatrixBase.hpp"

namespace AbstractLinAlgPack {

/** \brief Mix-in Interface for extracting the inverse cholesky factor of a dense symmetric
  * positive definite matrix.
  */
class MatrixExtractInvCholFactor
  : virtual public AbstractLinAlgPack::MatrixBase // doxygen needs full name
{
public:

  /** \brief Extract the inverse cholesly factor.
    *
    * Warning, the entire DMatrixSlice InvChol->gms() can be
    * used for workspace!
    */
  virtual void extract_inv_chol( DMatrixSliceTriEle* InvChol ) const = 0;	

};	// end class MatrixExtractInvCholFactor

}	// end namespace AbstractLinAlgPack

#endif // MATRIX_EXTRACT_INV_CHOL_FACTOR_H
