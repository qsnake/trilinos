/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// **********************************************************************/
 /* @HEADER@ */

#ifndef TSFINVERSEOPERATOR_IMPL_HPP
#define TSFINVERSEOPERATOR_IMPL_HPP

#include "SundanceDefs.hpp"
#include "SundanceTabs.hpp"
#include "TSFSolverState.hpp"
#include "TSFInverseOperatorDecl.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_VectorSpaceBase.hpp"

#include "Teuchos_RefCountPtr.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFSimpleTransposedOpImpl.hpp"
#endif


namespace TSFExtended
{
using Teuchos::RCP;

/*
 * Ctor with a linear operator and a solver specified.
 */
template <class Scalar> inline
InverseOperator<Scalar>::InverseOperator(const LinearOperator<Scalar>& op, 
  const LinearSolver<Scalar>& solver)
  : op_(op), solver_(solver) {;}


/* 
 * Compute alpha*M*x + beta*y, where M=*this.
 * @param M_trans specifies whether the operator is transposed:
 *                op(M) = M, for M_trans == NOTRANS
 *                op(M) = M', for M_trans == TRANS
 * @param x       vector of length this->domain()->dim()
 * @param y       vector of length this->range()->dim()
 * @param alpha   scalar multiplying M*x (default is 1.0)
 * @param beta    scalar multiplying y (default is 0.0)
 */
template <class Scalar> inline
void InverseOperator<Scalar>::generalApply(
  const Thyra::EOpTransp            M_trans
  ,const Thyra::VectorBase<Scalar>    &x
  ,Thyra::VectorBase<Scalar>          *y
  ,const Scalar            alpha 
  ,const Scalar            beta  
  ) const 
{
  
  Tabs tab(0);
  SUNDANCE_MSG2(this->verb(), tab << "InverseOperator::generalApply()");

  TEST_FOR_EXCEPTION(dynamic_cast<Thyra::ZeroLinearOpBase<Scalar>* >(op_.ptr().get()) != 0, std::runtime_error,
    "InverseOperator<Scalar>::apply() called on a ZeroOperator.");
  TEST_FOR_EXCEPTION(op_.domain().dim() != op_.range().dim(), std::runtime_error,
    "InverseOperator<Scalar>::apply() called on a non-square operator.");

  if (alpha==Teuchos::ScalarTraits<Scalar>::zero())
  {
    Ptr<VectorBase<Scalar> > yp(y);
    Vt_S(yp, beta);
  }
  else
  {
    Vector<Scalar> temp = createMember(*(x.space()));
    Vector<Scalar> result;
    assign(temp.ptr().ptr(), x);
    SolverState<Scalar> haveSoln;
    if (M_trans==Thyra::NOTRANS)
    {
      haveSoln = solver_.solve(op_, temp, result);
    }
    else
    {
      haveSoln = solver_.solve(op_.transpose(), temp, result);
    }
    TEST_FOR_EXCEPTION(haveSoln.finalState() != SolveConverged, 
      std::runtime_error,
      "InverseOperator<Scalar>::apply() " 
      << haveSoln.stateDescription());
    Vt_S(result.ptr().ptr(), alpha);
    Ptr<VectorBase<Scalar> > yp(y);
    V_StVpV(yp, beta, *y, *(result.ptr()));
  }      
  SUNDANCE_MSG2(this->verb(), tab << "done InverseOperator::generalApply()");
}


/** 
 * Return the domain of the operator. 
 */
template <class Scalar> inline
RCP<const Thyra::VectorSpaceBase<Scalar> > 
InverseOperator<Scalar>::domain() const 
{return op_.ptr()->domain();}
    

/** 
 * Return the range of the operator. 
 */
template <class Scalar> inline
RCP<const Thyra::VectorSpaceBase<Scalar> >
InverseOperator<Scalar>::range() const 
{return op_.ptr()->range();}


template <class Scalar> 
void InverseOperator<Scalar>::print(std::ostream& os) const
{
  Tabs tab(0);
  os << tab << "InverseOperator[" << std::endl;
  Tabs tab1;
  os << tab1 << "op=" << op_ << std::endl;
  os << tab << "]" << std::endl;
}


template <class Scalar> 
LinearOperator<Scalar> 
inverse(const LinearOperator<Scalar>& op, 
  const LinearSolver<Scalar>& solver)
{
  RCP<LinearOpBase<Scalar> > rtn 
    = rcp(new InverseOperator<Scalar>(op, solver));
  return rtn;
}

}

#endif
