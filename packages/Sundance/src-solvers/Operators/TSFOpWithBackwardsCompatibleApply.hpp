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

#ifndef TSFOPWITHBACKWARDSCOMPATIBLEAPPLY_HPP
#define TSFOPWITHBACKWARDSCOMPATIBLEAPPLY_HPP

#include "Thyra_VectorBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_AssertOp.hpp"
#include "Thyra_DefaultColumnwiseMultiVector.hpp"
#include "Teuchos_TestForException.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceNamedObject.hpp"
#include "SundanceOut.hpp"

namespace TSFExtended
{
using namespace Teuchos;
using namespace Thyra;
using namespace Sundance;

/** */
template <class Scalar>
class OpWithBackwardsCompatibleApply :
    public LinearOpBase<Scalar>,
    public DefaultObjectWithVerbosity
{
public:

  /** \brief . */
	bool opSupportedImpl(Thyra::EOpTransp M_trans) const
    {
      return (M_trans == Thyra::NOTRANS);
    }

  /** Thyra apply method */
  virtual void applyImpl(
    const Thyra::EOpTransp M_trans,
    const Thyra::MultiVectorBase<Scalar>    &X,
    const Ptr<Thyra::MultiVectorBase<Scalar> > &Y,
    const Scalar alpha,
    const Scalar beta
    ) const 
    {
      generalApply(M_trans, X, &*Y, alpha, beta);
    }

  /**
   * generalApply() applies either the operator or the transpose
   * according to the value of the transpose flag. This method is
   * backwards compatible with TSFCore-based code.
   */
  virtual void generalApply(const Thyra::EOpTransp M_trans,
    const Thyra::MultiVectorBase<Scalar>    &x,
    Thyra::MultiVectorBase<Scalar>          *y,
    const Scalar            alpha,
    const Scalar            beta) const 
    {

      using Teuchos::dyn_cast;

      THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
        "OpWithBackwardsCompatibleApply<Scalar>::applyImpl(...)",
        *this, M_trans, x, y );

/*
      const Thyra::DefaultColumnwiseMultiVector<Scalar> &xVec =
        dyn_cast<const Thyra::DefaultColumnwiseMultiVector<Scalar> >(x);
      
      Thyra::DefaultColumnwiseMultiVector<Scalar> &yVec = 
        dyn_cast<Thyra::DefaultColumnwiseMultiVector<Scalar> >(*y);
*/
      
      int nXCols = x.domain()->dim();
      
      for (int i=0; i < nXCols; i++)
      {
        generalApply(M_trans, *x.col(i), &*y->col(i), alpha, beta);
        //generalApply(M_trans, *xVec.col(i), &*yVec.col(i), alpha, beta);
      }
    }

/**
 *
 */
  virtual void generalApply(const Thyra::EOpTransp M_trans,
    const Thyra::VectorBase<Scalar>    &x,
    Thyra::VectorBase<Scalar>* y,
    const Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
    const Scalar beta  = Teuchos::ScalarTraits<Scalar>::zero()) const = 0 ;
};

}

#endif
