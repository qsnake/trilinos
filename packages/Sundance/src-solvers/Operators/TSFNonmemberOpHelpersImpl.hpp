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

#ifndef TSFNONMEMBEROPHELPERS_IMPL_HPP
#define TSFNONMEMBEROPHELPERS_IMPL_HPP


#include "TSFMultiVectorOperator.hpp"
#include "TSFCommonOperatorsDecl.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "SundanceOut.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFCommonOperatorsImpl.hpp"
#include "TSFLinearOperatorImpl.hpp"
#endif

namespace TSFExtended
{

template <class Scalar> inline
LinearOperator<Scalar> zeroOperator(
  const VectorSpace<Scalar>& domain,
  const VectorSpace<Scalar>& range)
{
  RCP<LinearOpBase<Scalar> > op 
    = rcp(new SimpleZeroOp<Scalar>(domain, range));

  return op;
}


template <class Scalar> inline
LinearOperator<Scalar> identityOperator(
  const VectorSpace<Scalar>& space)
{
  RCP<LinearOpBase<Scalar> > op 
    = rcp(new SimpleIdentityOp<Scalar>(space));

  return op;
}


template <class Scalar> inline
LinearOperator<Scalar> diagonalOperator(
  const Vector<Scalar>& vec)
{
  RCP<LinearOpBase<Scalar> > op 
    = rcp(new DefaultDiagonalLinearOp<Scalar>(vec.ptr()));

  return op;
}



template <class Scalar> inline
LinearOperator<Scalar> composedOperator(
  const Array<LinearOperator<Scalar> >& ops)
{
  /* We will strip out any identity operators, and if we find a zero
  * operator the whole works becomes a zero operator */ 
  Array<LinearOperator<Scalar> > strippedOps;

  for (int i=0; i<ops.size(); i++)
  {
    LinearOperator<Scalar> op_i = ops[i];

    /* if a factor is zero, the whole operator is
     * a zero operator */
    const SimpleZeroOp<Scalar>* zPtr 
      = dynamic_cast<const SimpleZeroOp<Scalar>*>(op_i.ptr().get());

    if (zPtr != 0) 
    {
      VectorSpace<Scalar> r = ops[0].range();
      VectorSpace<Scalar> d = ops[ops.size()-1].domain();
      return zeroOperator(d, r);
    }

    /* if a factor is the identity, skip it */
    const SimpleIdentityOp<Scalar>* IPtr 
      = dynamic_cast<const SimpleIdentityOp<Scalar>*>(op_i.ptr().get());  
    if (IPtr != 0) 
    {
      continue;
    }

    strippedOps.append(op_i);
  }
  
  TEST_FOR_EXCEPT(strippedOps.size() < 1);
  if (strippedOps.size()==1) return strippedOps[0];
  
  RCP<LinearOpBase<Scalar> > op 
    = rcp(new SimpleComposedOp<Scalar>(strippedOps));
  return op;
}





template <class Scalar> inline
LinearOperator<Scalar> addedOperator(
  const Array<LinearOperator<Scalar> >& ops)
{
  /* We will strip out any zero operators */
  Array<LinearOperator<Scalar> > strippedOps;

  for (int i=0; i<ops.size(); i++)
  {
    LinearOperator<Scalar> op_i = ops[i];

    /* Ignore any zero operators */
    const SimpleZeroOp<Scalar>* zPtr 
      = dynamic_cast<const SimpleZeroOp<Scalar>*>(op_i.ptr().get());

    if (zPtr != 0) continue;

    strippedOps.append(op_i);
  }
  
  TEST_FOR_EXCEPT(strippedOps.size() < 1);
  if (strippedOps.size()==1) return strippedOps[0];
  
  RCP<LinearOpBase<Scalar> > op 
    = rcp(new SimpleAddedOp<Scalar>(strippedOps));
  
  return op;
}



template <class Scalar> inline
LinearOperator<Scalar> scaledOperator(
  const Scalar& scale,
  const LinearOperator<Scalar>& op)
{
  RCP<LinearOpBase<Scalar> > A 
    = rcp(new DefaultScaledAdjointLinearOp<Scalar>(scale, Thyra::NOTRANS, op.ptr()));

  return A;
}



template <class Scalar> inline
LinearOperator<Scalar> scaledTransposedOperator(
  const Scalar& scale,
  const LinearOperator<Scalar>& op)
{
  RCP<LinearOpBase<Scalar> > A 
    = rcp(new DefaultScaledAdjointLinearOp<Scalar>(scale, Thyra::TRANS, op.ptr()));

  return A;
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

  
template <class Scalar> inline
LinearOperator<Scalar> multiVectorOperator(
  const Teuchos::Array<Vector<Scalar> >& cols,
  const VectorSpace<Scalar>& domain)
{
  RCP<LinearOpBase<Scalar> > A
    = rcp(new MultiVectorOperator<Scalar>(cols, domain));

  return A;
}

 
  
}


#endif
