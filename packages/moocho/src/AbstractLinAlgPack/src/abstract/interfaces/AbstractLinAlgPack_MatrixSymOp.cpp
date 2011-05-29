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

#include <assert.h>

#include "AbstractLinAlgPack_MatrixSymOp.hpp"
#include "AbstractLinAlgPack_EtaVector.hpp"

namespace AbstractLinAlgPack {

MatrixSymOp::mat_mswo_mut_ptr_t
MatrixSymOp::clone_mswo()
{
  return Teuchos::null;
}

MatrixSymOp::mat_mswo_ptr_t
MatrixSymOp::clone_mswo() const
{
  return Teuchos::null;
}

void MatrixSymOp::Mp_StPtMtP(
  MatrixSymOp* sym_lhs, value_type alpha
  , EMatRhsPlaceHolder dummy_place_holder
  , const GenPermMatrixSlice& gpms_rhs, BLAS_Cpp::Transp gpms_rhs_trans
  , value_type beta ) const
{
  TEST_FOR_EXCEPT(true); // ToDo: Implement!
}

void MatrixSymOp::Mp_StMtMtM(
  MatrixSymOp* sym_lhs, value_type alpha
  , EMatRhsPlaceHolder dummy_place_holder
  , const MatrixOp& mwo_rhs, BLAS_Cpp::Transp mwo_rhs_trans
  , value_type beta ) const
{
  TEST_FOR_EXCEPT(true); // ToDo: Implement!
}

// Overridden from MatrixOp


size_type MatrixSymOp::cols() const
{
  return this->rows();
}

const VectorSpace& MatrixSymOp::space_rows() const
{
  return this->space_cols();
}

MatrixSymOp::mat_mut_ptr_t
MatrixSymOp::clone()
{
  return clone_mswo();
}

MatrixSymOp::mat_ptr_t
MatrixSymOp::clone() const
{
  return clone_mswo();
}

}	// end namespace AbstractLinAlgPack 
