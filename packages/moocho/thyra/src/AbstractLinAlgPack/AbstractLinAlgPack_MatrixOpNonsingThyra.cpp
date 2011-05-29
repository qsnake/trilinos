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

#include "AbstractLinAlgPack_MatrixOpNonsingThyra.hpp"
#include "AbstractLinAlgPack_VectorMutableThyra.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace AbstractLinAlgPack {

// Constructors / Initializers

MatrixOpNonsingThyra::MatrixOpNonsingThyra()
{}

MatrixOpNonsingThyra::MatrixOpNonsingThyra(
  const Teuchos::RCP<const Thyra::LinearOpWithSolveBase<value_type> >  &thyra_linear_op_ns
  ,BLAS_Cpp::Transp                                                            thyra_linear_op_trans
  )
{
  this->initialize(thyra_linear_op_ns,thyra_linear_op_trans);
}

void MatrixOpNonsingThyra::initialize(
  const Teuchos::RCP<const Thyra::LinearOpWithSolveBase<value_type> >  &thyra_linear_op_ns
  ,BLAS_Cpp::Transp                                                            thyra_linear_op_trans
  )
{
  namespace mmp = MemMngPack;
  TEST_FOR_EXCEPTION(
    thyra_linear_op_ns.get()==NULL, std::invalid_argument
    ,"MatrixOpNonsingThyra::initialize(thyra_linear_op_ns): Error!"
    );
  MatrixOpThyra::initialize(thyra_linear_op_ns,thyra_linear_op_trans);
}

Teuchos::RCP<const Thyra::LinearOpWithSolveBase<value_type> > 
MatrixOpNonsingThyra::set_uninitialized()
{
  Teuchos::RCP<const Thyra::LinearOpWithSolveBase<value_type> >
    tmp_thyra_linear_op_ns = thyra_linear_op_ns();
  MatrixOpThyra::set_uninitialized();
  return tmp_thyra_linear_op_ns;
}

Teuchos::RCP<const Thyra::LinearOpWithSolveBase<value_type> >
MatrixOpNonsingThyra::thyra_linear_op_ns() const
{
  return Teuchos::rcp_dynamic_cast<const Thyra::LinearOpWithSolveBase<value_type> >(this->thyra_linear_op());
}

// Overridden from MatrixOp (needed to remove ambiguities)

MatrixOp::mat_mut_ptr_t
MatrixOpNonsingThyra::clone()
{
  return this->MatrixOpThyra::clone();
}

// Overridden from MatrixNonsing

void MatrixOpNonsingThyra::V_InvMtV(
  VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
  ,const Vector& v_rhs2
  ) const
{
  using Teuchos::dyn_cast;
  using BLAS_Cpp::trans_trans;
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    v_lhs==NULL, std::invalid_argument
    ,"MatrixOpThyra::Vp_StMtV(...): Error!"
    );
#endif
  *v_lhs = 0.0; // Must initialize before sending to solve(...)!
  VectorMutableThyra &v_thyra_lhs = dyn_cast<VectorMutableThyra>(*v_lhs);
  Teuchos::RCP<Thyra::VectorBase<value_type> > thyra_vec_lhs = v_thyra_lhs.set_uninitialized();
  Thyra::solve(
    *thyra_linear_op_ns()
    ,trans_trans(trans_rhs1,thyra_linear_op_trans())==BLAS_Cpp::no_trans ? Thyra::NOTRANS : Thyra::TRANS  // M_trans
    ,*dyn_cast<const VectorMutableThyra>(v_rhs2).thyra_vec()                                              // y
    ,thyra_vec_lhs.get()                                                                                  // x
    );
  v_thyra_lhs.initialize(thyra_vec_lhs);
}

void MatrixOpNonsingThyra::M_StInvMtM(
  MatrixOp* m_lhs, value_type alpha
  ,BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
  ) const
{
  MatrixNonsing::M_StInvMtM(m_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2); // ToDo: Specialize!
}

// Overridden from MatrixOpNonsing

MatrixOpNonsing::mat_mwons_ptr_t
MatrixOpNonsingThyra::clone_mwons() const
{
  return Teuchos::null; // ToDo: Add a clone function to Thyra::LinearOpWithSolveBase???
  //return Teuchos::rcp(new MatrixOpNonsingThyra(thyra_linear_op_ns()->clone_lows()));
}

} // end namespace AbstractLinAlgPack
