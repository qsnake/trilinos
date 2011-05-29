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

#include "AbstractLinAlgPack_MultiVectorMutable.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "RTOp_TOp_assign_scalar.h"
#include "RTOp_TOp_assign_vectors.h"
#include "RTOp_TOp_scale_vector.h"
#include "RTOpPack_RTOpC.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace {

// vector scalar assignment operator
RTOpPack::RTOpC& assign_scalar_op()
{
  static RTOpPack::RTOpC          assign_scalar_op_;
  return(assign_scalar_op_);
}
// vector assignment operator
static RTOpPack::RTOpC          assign_vec_op;
// scale vector
static RTOpPack::RTOpC          scale_vector_op;

// Simple class for an object that will initialize the operator objects
class init_rtop_server_t {
public:
  init_rtop_server_t() {
    // Vector scalar assignment operator
    TEST_FOR_EXCEPT(0!=RTOp_TOp_assign_scalar_construct(0.0,&assign_scalar_op().op()));
    // Vector assignment operator
    TEST_FOR_EXCEPT(0!=RTOp_TOp_assign_vectors_construct(&assign_vec_op.op()));
    // Operator scale_vector
    TEST_FOR_EXCEPT(0!=RTOp_TOp_scale_vector_construct(0.0,&scale_vector_op.op()));
  }
}; 

// When the program starts, this object will be created and the RTOp_Server object will
// be initialized before main() gets underway!
init_rtop_server_t  init_rtop_server;

} // end namespace

namespace AbstractLinAlgPack {

// Clone

MultiVectorMutable::multi_vec_mut_ptr_t
MultiVectorMutable::mv_clone()
{
  multi_vec_mut_ptr_t
    new_mv = this->space_cols().create_members(this->cols());
  const MultiVector*  multi_vecs[]      = { this };
  MultiVectorMutable* targ_multi_vecs[] = { new_mv.get() };
  AbstractLinAlgPack::apply_op(APPLY_BY_COL,assign_vec_op,1,multi_vecs,1,targ_multi_vecs,NULL);
  return new_mv;
}

// Sub-view methods

MultiVectorMutable::multi_vec_mut_ptr_t
MultiVectorMutable::mv_sub_view(const Range1D& row_rng, const Range1D& col_rng)
{
  TEST_FOR_EXCEPT(true); // ToDo: return a MultiVectorMutableSubView object.
  // Note that the MultiVectorMutableSubView class should derive from
  // MultiVectorSubView.
  return Teuchos::null;
}

// Overridden from MatrixOp

void MultiVectorMutable::zero_out()
{
  TEST_FOR_EXCEPT(0!=RTOp_TOp_assign_scalar_set_alpha(0.0,&assign_scalar_op().op()));
  MultiVectorMutable* targ_multi_vecs[] = { this };
  AbstractLinAlgPack::apply_op(APPLY_BY_COL,assign_scalar_op(),0,NULL,1,targ_multi_vecs,NULL);
}

void MultiVectorMutable::Mt_S( value_type alpha )
{
  if( alpha == 0.0 ) {
    TEST_FOR_EXCEPT(0!=RTOp_TOp_assign_scalar_set_alpha(alpha,&assign_scalar_op().op()));
    MultiVectorMutable* targ_multi_vecs[] = { this };
    AbstractLinAlgPack::apply_op(APPLY_BY_COL,assign_scalar_op(),0,NULL,1,targ_multi_vecs,NULL);
  }
  else if( alpha != 1.0 ) {
    TEST_FOR_EXCEPT(0!=RTOp_TOp_scale_vector_set_alpha(alpha,&scale_vector_op.op()));
    MultiVectorMutable* targ_multi_vecs[] = { this };
    AbstractLinAlgPack::apply_op(APPLY_BY_COL,scale_vector_op,0,NULL,1,targ_multi_vecs,NULL);
  }
}

MatrixOp& MultiVectorMutable::operator=(const MatrixOp& mwo_rhs)
{
  const MultiVector *mv_rhs = dynamic_cast<const MultiVector*>(&mwo_rhs);
  if(mv_rhs) {
    const MultiVector*  multi_vecs[]      = { mv_rhs };
    MultiVectorMutable* targ_multi_vecs[] = { this };
    AbstractLinAlgPack::apply_op(APPLY_BY_COL,assign_vec_op,1,multi_vecs,1,targ_multi_vecs,NULL);
  }
  else {
    TEST_FOR_EXCEPT(true); // ToDo: Get column by column or row by row
  }
  return *this;
}

MatrixOp::mat_mut_ptr_t
MultiVectorMutable::clone()
{
  return this->mv_clone();
}

bool MultiVectorMutable::Mp_StM(
  MatrixOp* mwo_lhs, value_type alpha
  ,BLAS_Cpp::Transp trans_rhs
  ) const
{
  return false; // ToDo: Specialize!
}

bool MultiVectorMutable::Mp_StM(
  value_type alpha,const MatrixOp& M_rhs, BLAS_Cpp::Transp trans_rhs
  )
{
  return false; // ToDo: Specialize!
}

// Overridden form MultiVector

MultiVector::multi_vec_ptr_t MultiVectorMutable::mv_clone() const
{
  return const_cast<MultiVectorMutable*>(this)->mv_clone();
}

MultiVectorMutable::vec_ptr_t MultiVectorMutable::col(index_type j) const
{
  return const_cast<MultiVectorMutable*>(this)->col(j);
}

MultiVectorMutable::vec_ptr_t MultiVectorMutable::row(index_type i) const
{
  return const_cast<MultiVectorMutable*>(this)->row(i);
}

MultiVectorMutable::vec_ptr_t MultiVectorMutable::diag(int k) const
{
  return const_cast<MultiVectorMutable*>(this)->diag(k);
}

MultiVectorMutable::multi_vec_ptr_t
MultiVectorMutable::mv_sub_view(const Range1D& row_rng, const Range1D& col_rng) const
{
  return const_cast<MultiVectorMutable*>(this)->mv_sub_view(row_rng,col_rng);
}

} // end namespace AbstractLinAlgPack
