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

#ifndef TSF_SIMPLE_ADDED_OP_IMPL_HPP
#define TSF_SIMPLE_ADDED_OP_IMPL_HPP



#include "TSFSimpleAddedOpDecl.hpp"
#include "TSFSimpleZeroOpDecl.hpp"
#include "TSFLinearCombinationImpl.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "Teuchos_Array.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFSimpleZeroOpImpl.hpp"
#include "TSFLinearOperatorImpl.hpp"
#include "TSFSimplifiedLinearOpBaseImpl.hpp"
#include "TSFSimpleTransposedOpImpl.hpp"
#endif


namespace TSFExtended
{
using namespace Teuchos;
using namespace Sundance;


/*
 * Represent a sum of operators A_0 + A_1 + ... + A_n.
 */
template <class Scalar> inline
SimpleAddedOp<Scalar>::SimpleAddedOp(
  const Array<LinearOperator<Scalar> >& ops)
  : SimplifiedLinearOpWithSpaces<Scalar>(
    ops[0].domain(), ops[0].range()
    ) 
  , ops_(ops)
{
  TEST_FOR_EXCEPT(ops_.size() <= 1);
  for (int i=1; i<ops_.size(); i++)
  {
    TEST_FOR_EXCEPT(!(ops[i].range() == ops[0].range()));
    TEST_FOR_EXCEPT(!(ops[i].domain() == ops[0].domain()));
  }
}
  
/* */
template <class Scalar> inline
void SimpleAddedOp<Scalar>::applyOp(const Thyra::EOpTransp M_trans,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab(0);
  SUNDANCE_MSG2(this->verb(), tab << "SimpleAddedOp::applyOp()");

  Vector<Scalar> tmp=out.copy();
  tmp.zero();
  for (int i=0; i<ops_.size(); i++)
  {
    Tabs tab1;
    Out::os() << tab1 << "applying term i=" << i << " of " 
              << ops_.size() << std::endl;
    if (M_trans == Thyra::NOTRANS)
      tmp = tmp + ops_[i] * in;
    else if (M_trans == Thyra::TRANS)
      tmp = tmp + ops_[i].transpose() * in;
    else 
      TEST_FOR_EXCEPT(M_trans != Thyra::TRANS && M_trans != Thyra::NOTRANS);
  }
  out.acceptCopyOf(tmp);

  SUNDANCE_MSG2(this->verb(), tab << "done SimpleAddedOp::applyOp()");
}
  
/* */
template <class Scalar> inline
std::string SimpleAddedOp<Scalar>::description() const 
{
  std::string rtn="(";
  for (int i=0; i<ops_.size(); i++)
  {
    if (i > 0) rtn += "+";
    rtn += ops_[i].description();
  }
  rtn += ")";
  return rtn;
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
LinearOperator<Scalar> operator+(const LinearOperator<Scalar>& A,
  const LinearOperator<Scalar>& B)
{
  return addedOperator(Array<LinearOperator<Scalar> >(tuple(A, B)));
}

}

#endif
