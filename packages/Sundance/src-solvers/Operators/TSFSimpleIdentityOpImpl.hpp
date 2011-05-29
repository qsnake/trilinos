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

#ifndef TSF_SIMPLE_IDENTITY_OP_IMPL_HPP
#define TSF_SIMPLE_IDENTITY_OP_IMPL_HPP



#include "TSFSimpleIdentityOpDecl.hpp"
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


/* ---- Identity op ------- */

template <class Scalar> inline
SimpleIdentityOp<Scalar>::SimpleIdentityOp(const VectorSpace<Scalar>& space)
  : SimplifiedLinearOpWithSpaces<Scalar>(space, space) {}


template <class Scalar> inline
void SimpleIdentityOp<Scalar>::applyOp(const Thyra::EOpTransp M_trans,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab(0);
  SUNDANCE_MSG2(this->verb(), tab << "SimpleZeroOp::applyOp()");
  out.acceptCopyOf(in);
  SUNDANCE_MSG2(this->verb(), tab << "done SimpleIdentityOp::applyOp()");
}

template <class Scalar> inline
std::string SimpleIdentityOp<Scalar>::description() const 
{return "I(" + this->domain()->description() + ")";}



template <class Scalar> inline
LinearOperator<Scalar> identityOperator(
  const VectorSpace<Scalar>& space)
{
  RCP<LinearOpBase<Scalar> > op 
    = rcp(new SimpleIdentityOp<Scalar>(space));

  return op;
}



}

#endif
