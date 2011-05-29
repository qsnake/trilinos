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

#ifndef TSF_SIMPLE_COMPOSED_OP_IMPL_HPP
#define TSF_SIMPLE_COMPOSED_OP_IMPL_HPP



#include "TSFSimpleComposedOpDecl.hpp"
#include "TSFSimpleIdentityOpDecl.hpp"
#include "TSFSimpleZeroOpDecl.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#include "TSFSimpleIdentityOpImpl.hpp"
#include "TSFSimpleZeroOpImpl.hpp"
#include "TSFSimplifiedLinearOpBaseImpl.hpp"
#endif


namespace TSFExtended
{
using namespace Teuchos;
using namespace Sundance;



/*
 * ------ composed operator  
 */
template <class Scalar> inline
SimpleComposedOp<Scalar>::SimpleComposedOp(const Array<LinearOperator<Scalar> >& ops)
  : SimplifiedLinearOpWithSpaces<Scalar>(
    ops[ops.size()-1].domain(), ops[0].range()
    ) 
  , ops_(ops)
{
  TEST_FOR_EXCEPT(ops_.size() <= 1);
  for (int i=1; i<ops_.size(); i++)
  {
    TEST_FOR_EXCEPT(!(ops[i].range() == ops[i-1].domain()));
  }
}
  


template <class Scalar> inline
void SimpleComposedOp<Scalar>::applyOp(const Thyra::EOpTransp M_trans,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab(0);
  SUNDANCE_MSG2(this->verb(), tab << "SimpleComposedOp::applyOp()");
  if (M_trans == Thyra::NOTRANS)
  {
    Vector<Scalar> tmpIn = in.copy();
    for (int i=0; i<ops_.size(); i++)
    {
      Tabs tab1;
      Vector<Scalar> tmpOut;
      int j = ops_.size()-1-i;
      SUNDANCE_MSG2(this->verb(), tab1 << "applying op #" << j 
        << " of " << ops_.size());
      ops_[j].apply(tmpIn, tmpOut);
      tmpIn = tmpOut;
    }
    out.acceptCopyOf(tmpIn);
  }

  else if (M_trans == Thyra::TRANS)
  {
    Vector<Scalar> tmpIn = in.copy();
    for (int i=0; i<ops_.size(); i++)
    {
      Tabs tab1;
      Vector<Scalar> tmpOut;
      SUNDANCE_MSG2(this->verb(), tab1 << "applying transpose op #" << i
        << " of " << ops_.size());
      ops_[i].applyTranspose(tmpIn, tmpOut);
      tmpIn = tmpOut;
    }
    out.acceptCopyOf(tmpIn);
  }
  else
  {
    TEST_FOR_EXCEPT(M_trans != Thyra::TRANS && M_trans != Thyra::NOTRANS);
  }
  SUNDANCE_MSG2(this->verb(), tab << "done SimpleComposedOp::applyOp()");
}
  



template <class Scalar> inline
std::string SimpleComposedOp<Scalar>::description() const 
{
  std::string rtn="(";
  for (int i=0; i<ops_.size(); i++)
  {
    if (i > 0) rtn += "*";
    rtn += ops_[i].description();
  }
  rtn += ")";
  return rtn;
}


template <class Scalar> inline
void SimpleComposedOp<Scalar>::print(std::ostream& os) const 
{
  Tabs tab(0);
  os << tab << "ComposedOperator[" << std::endl;
  for (int i=0; i<ops_.size(); i++)
  {
    Tabs tab1;
    os << tab1 << "factor #" << i << std::endl;
    Tabs tab2;
    os << tab2 << ops_[i].description() << std::endl;
    os << std::endl;
  }
  os << tab << "]" <<  std::endl;
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
LinearOperator<Scalar> operator*(const LinearOperator<Scalar>& A, 
  const LinearOperator<Scalar>& B)
{
  return composedOperator(Array<LinearOperator<Scalar> >(tuple(A,B)));
}
  

}

#endif
