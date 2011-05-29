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

#ifndef TSF_SIMPLE_SCALED_OP_IMPL_HPP
#define TSF_SIMPLE_SCALED_OP_IMPL_HPP



#include "TSFSimpleScaledOpDecl.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#include "TSFSimplifiedLinearOpBaseImpl.hpp"
#endif


namespace TSFExtended
{
using namespace Teuchos;
using namespace Sundance;




/*
 * --- scaled op
 */

template <class Scalar> inline
SimpleScaledOp<Scalar>::SimpleScaledOp(const Scalar& alpha,
  const LinearOperator<Scalar>& A)
  : SimplifiedLinearOpWithSpaces<Scalar>(
    A.domain(), A.range()
    ) 
  , alpha_(alpha), A_(A)
{}
  
/* */
template <class Scalar> inline
void SimpleScaledOp<Scalar>::applyOp(const Thyra::EOpTransp M_trans,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab(0);
  SUNDANCE_MSG2(this->verb(), tab << "SimpleScaledOp::applyOp()");

  if (M_trans == Thyra::NOTRANS)
    A_.apply(in, out);
  else if (M_trans == Thyra::TRANS)
    A_.applyTranspose(in, out);
  else 
    TEST_FOR_EXCEPT(M_trans !=Thyra::TRANS && M_trans != Thyra::NOTRANS);

  out.scale(alpha_);

  SUNDANCE_MSG2(this->verb(), tab << "done SimpleScaledOp::applyOp()");
}
  
/* */
template <class Scalar> inline
std::string SimpleScaledOp<Scalar>::description() const 
{
  return "ScaledOp[alpha="  + Teuchos::toString(alpha_)
    + ", " + A_.description() + "]";
}


/* */
template <class Scalar> inline
void SimpleScaledOp<Scalar>::print(std::ostream& os) const 
{
  Tabs tab(0);
  os << tab << "ScaledOp[" << std::endl;
  Tabs tab1;
  os << tab1 << "scale = " << alpha_ << std::endl;
  os << tab1 << "operator = " << A_.description() << std::endl;
  os << tab << "]" << std::endl;
}



template <class Scalar> inline
LinearOperator<Scalar> scaledOperator(
  const Scalar& scale,
  const LinearOperator<Scalar>& op)
{
  RCP<LinearOpBase<Scalar> > A 
    = rcp(new SimpleScaledOp<Scalar>(scale, op));

  return A;
}


template <class Scalar> inline
LinearOperator<Scalar> operator*(const Scalar& a, const LinearOperator<Scalar>& A)
{
  return scaledOperator(a, A);
}
  

}

#endif
