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

#include "AbstractLinAlgPack_MatrixOpThyra.hpp"
#include "AbstractLinAlgPack_VectorMutableThyra.hpp"
#include "AbstractLinAlgPack_ThyraAccessors.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace AbstractLinAlgPack {

// Constructors / Initializers

MatrixOpThyra::MatrixOpThyra()
{}

MatrixOpThyra::MatrixOpThyra(
  const Teuchos::RCP<const Thyra::LinearOpBase<value_type> >   &thyra_linear_op
  ,BLAS_Cpp::Transp                                                    thyra_linear_op_trans
  )
{
  this->initialize(thyra_linear_op,thyra_linear_op_trans);
}

void MatrixOpThyra::initialize(
  const Teuchos::RCP<const Thyra::LinearOpBase<value_type> >   &thyra_linear_op
  ,BLAS_Cpp::Transp                                                    thyra_linear_op_trans
  )
{
  namespace mmp = MemMngPack;
  TEST_FOR_EXCEPTION(
    thyra_linear_op.get()==NULL, std::invalid_argument
    ,"MatrixOpThyra::initialize(thyra_linear_op): Error!"
    );
  const bool adjointSupported = ( ::Thyra::opSupported(*thyra_linear_op,Thyra::NOTRANS) && ::Thyra::opSupported(*thyra_linear_op,Thyra::TRANS) );
  TEST_FOR_EXCEPTION(
    !adjointSupported, std::invalid_argument
    ,"MatrixOpThyra::initialize(thyra_linear_op): Error, the operator opSupported(thyra_linear_op,transp) must return true "
    "for both values of transp==NOTRANS and transp=TRANS!"
    );
  thyra_linear_op_       = thyra_linear_op;
  thyra_linear_op_trans_ = thyra_linear_op_trans;
  space_cols_.initialize(thyra_linear_op_->range());
  space_rows_.initialize(thyra_linear_op_->domain());
}

Teuchos::RCP<const Thyra::LinearOpBase<value_type> > 
MatrixOpThyra::set_uninitialized()
{
  Teuchos::RCP<const Thyra::LinearOpBase<value_type> > tmp_thyra_linear_op = thyra_linear_op_;
  thyra_linear_op_       = Teuchos::null;
  thyra_linear_op_trans_ = BLAS_Cpp::no_trans;
  space_cols_.set_uninitialized();
  space_rows_.set_uninitialized();
  return tmp_thyra_linear_op;
}

// Overridden from MatrixBase

const VectorSpace&
MatrixOpThyra::space_cols() const
{
  return space_cols_;
}

const VectorSpace&
MatrixOpThyra::space_rows() const
{
  return space_rows_;
}

// Overridden from MatrixOp

MatrixOp::mat_mut_ptr_t
MatrixOpThyra::clone()
{
  return Teuchos::rcp(new MatrixOpThyra(thyra_linear_op_->clone()));
}

MatrixOp& MatrixOpThyra::operator=(const MatrixOp& mwo_rhs)
{
  using Teuchos::dyn_cast;
  const MatrixOpThyra &mwo_rhs_thyra = dyn_cast<const MatrixOpThyra>(mwo_rhs);
  thyra_linear_op_       = mwo_rhs_thyra.thyra_linear_op_;  // ToDo: Clone this!
  thyra_linear_op_trans_ = mwo_rhs_thyra.thyra_linear_op_trans_;
  space_cols_           = mwo_rhs_thyra.space_cols_;  // ToDo: Clone this!
  space_rows_           = mwo_rhs_thyra.space_rows_;  // ToDo: Clone this!
  return *this;
}

void MatrixOpThyra::Vp_StMtV(
  VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
  ,const Vector& v_rhs2, value_type beta
  ) const
{
  using BLAS_Cpp::trans_trans;
  // Get Thyra views of the vectors
  Teuchos::RCP<const Thyra::VectorBase<value_type> > thyra_vec_rhs2;
  get_thyra_vector( BLAS_Cpp::no_trans==trans_rhs1 ? space_rows_ : space_cols_, v_rhs2, &thyra_vec_rhs2 );
  Teuchos::RCP<Thyra::VectorBase<value_type> > thyra_vec_lhs;
  get_thyra_vector( BLAS_Cpp::no_trans==trans_rhs1 ? space_cols_ : space_rows_, v_lhs, &thyra_vec_lhs );
  // Perform the multiplication
  ::Thyra::apply(
    *thyra_linear_op_
    ,trans_trans(trans_rhs1,thyra_linear_op_trans())==BLAS_Cpp::no_trans ? Thyra::NOTRANS : Thyra::TRANS  // M_trans
    ,*thyra_vec_rhs2                                                                                      // x
    ,&*thyra_vec_lhs                                                                                      // y
    ,alpha                                                                                                // alpha
    ,beta                                                                                                // beta
    );
  // Free/commit Thyra vector views
  free_thyra_vector( BLAS_Cpp::no_trans==trans_rhs1 ? space_rows_ : space_cols_, v_rhs2, &thyra_vec_rhs2 );
  commit_thyra_vector( BLAS_Cpp::no_trans==trans_rhs1 ? space_cols_ : space_rows_, v_lhs, &thyra_vec_lhs );
}

bool MatrixOpThyra::Mp_StMtM(
  MatrixOp* mwo_lhs, value_type alpha
  ,BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
  ,value_type beta
  ) const
{
  return false; // ToDo: Specialize!
}

} // end namespace AbstractLinAlgPack
