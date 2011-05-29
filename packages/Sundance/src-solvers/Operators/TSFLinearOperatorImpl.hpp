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

#ifndef TSFLINEAROPERATORIMPL_HPP
#define TSFLINEAROPERATORIMPL_HPP

#include "SundanceDefs.hpp"
#include "TSFLinearOperatorDecl.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "TSFVectorDecl.hpp"
#include "TSFVectorSpaceDecl.hpp"
#include "TSFInverseOperatorDecl.hpp"
#include "TSFSimpleTransposedOpDecl.hpp"
#include "TSFBlockOperatorBaseDecl.hpp"
#include "TSFVectorType.hpp"
#include "SundanceOut.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorImpl.hpp"
#endif



using namespace TSFExtended;
using namespace Teuchos;
using namespace Sundance;

template <class Scalar>
class InverseOperator;


//=======================================================================
template <class Scalar>
LinearOperator<Scalar>::LinearOperator() 
  : Handle<Thyra::LinearOpBase<Scalar> >(), verb_(0) {;}


//=======================================================================
template <class Scalar>
LinearOperator<Scalar>::LinearOperator(const RCP<Thyra::LinearOpBase<Scalar> >& smartPtr) 
  : Handle<Thyra::LinearOpBase<Scalar> >(smartPtr), verb_(0) {;}




//=======================================================================
template <class Scalar> inline 
void LinearOperator<Scalar>::apply(const Vector<Scalar>& in,
  Vector<Scalar>& out,
  const Scalar& alpha,
  const Scalar& beta) const
{
  Tabs tab(0);
  SUNDANCE_MSG1(this->verb(), tab << "Operator=" << this->description()
    << ",  calling apply() function");
  Tabs tab1;
  SUNDANCE_MSG1(this->verb(), tab1 << "alpha=" << alpha);
  SUNDANCE_MSG1(this->verb(), tab1 << "beta=" << beta);
  
  if (this->verb() > 2)
  {
    Tabs tab2;
    Out::os() << tab2 << "input vector = " << in << std::endl;
  }
  else if (this->verb() > 1)
  {
    Tabs tab2;
    Out::os() << tab2 << "input vector = " << in.description() << std::endl;
  }

  /* the result vector might not be initialized. If it's null,
   * create a new vector in the range space */
  if (out.ptr().get()==0)
  {
    Tabs tab2;
    SUNDANCE_MSG3(this->verb(), tab2 << "allocating output vector");
    out = this->range().createMember();
  }
  else
  {
    Tabs tab2;
    SUNDANCE_MSG3(this->verb(), tab2 << "using preallocated output vector");
  }

  this->ptr()->apply(Thyra::NOTRANS, *(in.ptr().ptr()),
    out.ptr().ptr(), alpha, beta);

  
  
  if (this->verb() > 2)
  {
    Tabs tab2;
    Out::os() << tab2 << "output vector = " << out << std::endl;
  }
  else if (this->verb() > 1)
  {
    Tabs tab2;
    Out::os() << tab2 << "output vector = " << out.description() << std::endl;
  }

  SUNDANCE_MSG1(this->verb(), tab << "Operator=" << this->description()
    << ",  done with apply() function");
  
}




//=======================================================================
template <class Scalar> inline 
void LinearOperator<Scalar>::applyTranspose(const Vector<Scalar>& in,
  Vector<Scalar>& out,
  const Scalar& alpha,
  const Scalar& beta) const
{
  Tabs tab(0);
  SUNDANCE_MSG1(this->verb(), tab << "Operator=" << this->description()
    << ",  calling applyTranspose() function");
  Tabs tab1;
  SUNDANCE_MSG1(this->verb(), tab1 << "alpha=" << alpha);
  SUNDANCE_MSG1(this->verb(), tab1 << "beta=" << beta);
  
  if (this->verb() > 2)
  {
    Tabs tab2;
    Out::os() << tab2 << "input vector = " << in << std::endl;
  }
  else if (this->verb() > 1)
  {
    Tabs tab2;
    Out::os() << tab2 << "input vector = " << in.description() << std::endl;
  }


  /* the result vector might not be initialized. If it's null,
   * create a new vector in the range space */
  if (out.ptr().get()==0)
  {
    Tabs tab2;
    SUNDANCE_MSG3(this->verb(), tab2 << "allocating output vector");
    out = this->domain().createMember();
  }
  else
  {
    Tabs tab2;
    SUNDANCE_MSG3(this->verb(), tab2 << "using preallocated output vector");
  }

  this->ptr()->apply(Thyra::TRANS, *(in.ptr().ptr()),
    out.ptr().ptr(), alpha, beta);

  
  
  if (this->verb() > 2)
  {
    Tabs tab2;
    Out::os() << tab2 << "output vector = " << out << std::endl;
  }
  else if (this->verb() > 1)
  {
    Tabs tab2;
    Out::os() << tab2 << "output vector = " << out.description() << std::endl;
  }

  SUNDANCE_MSG1(this->verb(), tab << "Operator=" << this->description()
    << ",  done with applyTranpose() function");
  
}


//=======================================================================
template <class Scalar>
RCP<Time>& LinearOperator<Scalar>::opTimer()
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("Low-level vector operations");
  return rtn;
}

//=======================================================================
template <class Scalar>
LinearOperator<Scalar> LinearOperator<Scalar>::transpose() const
{
  LinearOperator<Scalar> op = transposedOperator(*this);
  return op;
}





//=======================================================================
template <class Scalar>
RCP<LoadableMatrix<Scalar> > LinearOperator<Scalar>::matrix()
{
  RCP<LoadableMatrix<Scalar> > rtn 
    = rcp_dynamic_cast<LoadableMatrix<Scalar> >(this->ptr());
  return rtn;
}

//=======================================================================
template <class Scalar>
void LinearOperator<Scalar>::getRow(const int& row, 
  Teuchos::Array<int>& indices, 
  Teuchos::Array<Scalar>& values) const
{
  const RowAccessibleOp<Scalar>* val = 
    dynamic_cast<const RowAccessibleOp<Scalar>* >(this->ptr().get());
  TEST_FOR_EXCEPTION(val == 0, std::runtime_error, 
    "Operator not row accessible; getRow() not defined.");
  val->getRow(row, indices, values);
}

//=============================================================================
template <class Scalar>
int LinearOperator<Scalar>::numBlockRows() const
{
  const BlockOperatorBase<Scalar>* b = dynamic_cast<const BlockOperatorBase<Scalar>* >(this->ptr().get());
  if (b==0) return 1;
  return b->numBlockRows(); 
}

//=============================================================================
template <class Scalar>
int LinearOperator<Scalar>::numBlockCols() const
{
  const BlockOperatorBase<Scalar>* b = dynamic_cast<const BlockOperatorBase<Scalar>* >(this->ptr().get());
  if (b==0) return 1;
  return b->numBlockCols(); 
}


//=============================================================================
template <class Scalar>
const VectorSpace<Scalar> 
LinearOperator<Scalar>::range() const
{return this->ptr()->range();}
  

//=============================================================================
template <class Scalar>
void LinearOperator<Scalar>::setBlock(int i, int j, 
  const LinearOperator<Scalar>& sub) 
{
  SetableBlockOperatorBase<Scalar>* b = 
    dynamic_cast<SetableBlockOperatorBase<Scalar>* >(this->ptr().get());
  
  TEST_FOR_EXCEPTION(b == 0, std::runtime_error, 
    "Can't call setBlock since operator not SetableBlockOperatorBase");

  b->setBlock(i, j, sub);
} 



//=============================================================================
template <class Scalar>
const  VectorSpace<Scalar> 
LinearOperator<Scalar>::domain() const 
{return this->ptr()->domain();}



//=============================================================================
template <class Scalar>
LinearOperator<Scalar> LinearOperator<Scalar>::getBlock(const int &i, 
  const int &j) const 
{
  const BlockOperatorBase<Scalar>* b = 
    dynamic_cast<const BlockOperatorBase<Scalar>* >(this->ptr().get());
  
  if (b==0)
  {
    TEST_FOR_EXCEPTION(i != 0 || j != 0, std::runtime_error, 
      "nonzero block index (" << i << "," << j << ") into "
      "non-block operator");
    return *this;
  }
  return b->getBlock(i, j);
}


//=============================================================================
template <class Scalar>
LinearOperator<Scalar> LinearOperator<Scalar>::getNonconstBlock(const int &i, 
  const int &j) 
{
  BlockOperatorBase<Scalar>* b = 
    dynamic_cast<BlockOperatorBase<Scalar>* >(this->ptr().get());
  
  if (b==0)
  {
    TEST_FOR_EXCEPTION(i != 0 || j != 0, std::runtime_error, 
      "nonzero block index (" << i << "," << j << ") into "
      "non-block operator");
    return *this;
  }
  return b->getNonconstBlock(i, j);
}

 

//=============================================================================
template <class Scalar>
void LinearOperator<Scalar>::endBlockFill() 
{
  SetableBlockOperatorBase<Scalar>* b = 
    dynamic_cast<SetableBlockOperatorBase<Scalar>* >(this->ptr().get());
  
  TEST_FOR_EXCEPTION(b == 0, std::runtime_error, 
    "Can't call endBlockFill because operator is not a SetableBlockOperator");

  
  b->endBlockFill();
} 






#endif
