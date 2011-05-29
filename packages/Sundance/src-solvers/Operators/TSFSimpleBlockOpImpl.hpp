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

#ifndef TSF_SIMPLEBLOCKOP_IMPL_HPP
#define TSF_SIMPLEBLOCKOP_IMPL_HPP

#include "SundanceTabs.hpp"
#include "SundanceExceptions.hpp"
#include "TSFSimpleBlockOpDecl.hpp"
#include "TSFSimpleZeroOpDecl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#include "TSFSimpleZeroOpImpl.hpp"
#include "TSFSimplifiedLinearOpBaseImpl.hpp"
#endif




namespace TSFExtended
{
using namespace Teuchos;
using namespace Sundance;


/* ---- Simplified linear op with spaces ------- */

template <class Scalar> inline
SimpleBlockOp<Scalar>::SimpleBlockOp(
  const VectorSpace<Scalar>& domain,
  const VectorSpace<Scalar>& range)
  : SimplifiedLinearOpWithSpaces<Scalar>(domain, range), blocks_(range.numBlocks())
{
  for (int i=0; i<blocks_.size(); i++) 
  {
    blocks_[i] = Array<LinearOperator<Scalar> >(domain.numBlocks());
    for (int j=0; j<blocks_[i].size(); j++)
    {
      blocks_[i][j] = zeroOperator(domain.getBlock(j), range.getBlock(i));
    }
  }
}

template <class Scalar> inline
int SimpleBlockOp<Scalar>::numBlockRows() const
{
  return blocks_.size();
}

template <class Scalar> inline
int SimpleBlockOp<Scalar>::numBlockCols() const
{
  return blocks_[0].size();
}

template <class Scalar> inline
const LinearOperator<Scalar>& SimpleBlockOp<Scalar>::getBlock(int i, int j) const 
{
  return blocks_[i][j];
}

template <class Scalar> inline
LinearOperator<Scalar> SimpleBlockOp<Scalar>::getNonconstBlock(int i, int j) 
{
  return blocks_[i][j];
}


template <class Scalar> inline
void SimpleBlockOp<Scalar>::setBlock(int i, int j, 
  const LinearOperator<Scalar>& Aij) 
{
  blocks_[i][j] = Aij;
}

template <class Scalar> inline
void SimpleBlockOp<Scalar>::applyOp(const Thyra::EOpTransp M_trans,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab(0);
  SUNDANCE_MSG2(this->verb(), tab << "SimpleBlockOp::applyOp()");
  if (M_trans == Thyra::NOTRANS)
  {
    out.zero();
    for (int i=0; i<this->numBlockRows(); i++)
    {
      for (int j=0; j<this->numBlockCols(); j++)
      {
        Vector<Scalar> tmp; 
        blocks_[i][j].apply(in.getBlock(j), tmp);
        out.getBlock(i).update(1.0, tmp);
      }
    }
  }

  else if (M_trans == Thyra::TRANS)
  {
    for (int i=0; i<this->numBlockCols(); i++)
    {
      out.zero();
      for (int j=0; j<this->numBlockRows(); j++)
      {
        Vector<Scalar> tmp;
        blocks_[j][i].applyTranspose(in.getBlock(j),tmp);
        out.getBlock(i).update(1.0, tmp);
      }
    }
  }
  else
  {
    TEST_FOR_EXCEPT(M_trans != Thyra::TRANS && M_trans != Thyra::NOTRANS);
  }
  SUNDANCE_MSG2(this->verb(), tab << "done SimpleBlockOp::applyOp()");
}


template <class Scalar> inline
LinearOperator<Scalar> makeBlockOperator(
  const VectorSpace<Scalar>& domain,
  const VectorSpace<Scalar>& range
  )
{
  RCP<SimpleBlockOp<Scalar> > b = 
    rcp(new SimpleBlockOp<Scalar>(domain, range));
  RCP<Thyra::LinearOpBase<Scalar> > p = b;
  return p;
}



}

#endif
