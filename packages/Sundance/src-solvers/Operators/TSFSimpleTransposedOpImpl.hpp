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

#ifndef TSF_SIMPLE_TRANSPOSED_OP_IMPL_HPP
#define TSF_SIMPLE_TRANSPOSED_OP_IMPL_HPP



#include "TSFSimpleTransposedOpDecl.hpp"
#include "TSFSimpleZeroOpDecl.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#include "TSFSimpleZeroOpImpl.hpp"
#include "TSFSimplifiedLinearOpBaseImpl.hpp"
#endif


namespace TSFExtended
{
using namespace Teuchos;
using namespace Sundance;



/*
 * --- transposed op
 */

template <class Scalar> inline
SimpleTransposedOp<Scalar>::SimpleTransposedOp(const LinearOperator<Scalar>& A)
  : SimplifiedLinearOpWithSpaces<Scalar>(
    A.range(), A.domain()
    ) 
  , A_(A)
{}
  
/* */
template <class Scalar> inline
void SimpleTransposedOp<Scalar>::applyOp(const Thyra::EOpTransp M_trans,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab(0);
  SUNDANCE_MSG2(this->verb(), tab << "SimpleTransposedOp::applyOp()");

  if (M_trans == Thyra::NOTRANS)
    A_.applyTranspose(in, out);
  else if (M_trans == Thyra::TRANS)
    A_.apply(in, out);
  else 
    TEST_FOR_EXCEPT(M_trans !=Thyra::TRANS && M_trans != Thyra::NOTRANS);

  SUNDANCE_MSG2(this->verb(), tab << "done SimpleTransposedOp::applyOp()");
}
  
/* */
template <class Scalar> inline
std::string SimpleTransposedOp<Scalar>::description() const 
{
  return "(" + A_.description() + "^T)";
}



template <class Scalar> inline
LinearOperator<Scalar> transposedOperator(
  const LinearOperator<Scalar>& op)
{

  /* If the operator is a transpose, return the untransposed op */
  const SimpleTransposedOp<Scalar>* tPtr
    = dynamic_cast<const SimpleTransposedOp<Scalar>*>(op.ptr().get());
  if (tPtr)
  {
    return tPtr->op();
  }

  /* If the operator is zero, return a transposed zero */
  const SimpleZeroOp<Scalar>* zPtr 
    = dynamic_cast<const SimpleZeroOp<Scalar>*>(op.ptr().get());

  if (zPtr != 0) 
  {
    VectorSpace<Scalar> r = op.range();
    VectorSpace<Scalar> d = op.domain();
    return zeroOperator(r, d);
  }


  /* Return a transposed operator */
  RCP<LinearOpBase<Scalar> > A
    = rcp(new SimpleTransposedOp<Scalar>(op));
      
  return A;
}

  


}

#endif
