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

#ifndef TSFSIMPLIFIEDOPBASE_IMPL_HPP
#define TSFSIMPLIFIEDOPBASE_IMPL_HPP

#include "SundanceDefs.hpp"
#include "SundanceOut.hpp"
#include "TSFSimplifiedLinearOpBaseDecl.hpp"
#include "TSFVectorDecl.hpp"




namespace TSFExtended
{
using namespace Teuchos;
using namespace Thyra;

template <class Scalar>
void SimplifiedLinearOpBase<Scalar>
::generalApply(const Thyra::EOpTransp M_trans,
  const Thyra::VectorBase<Scalar>    &x,
  Thyra::VectorBase<Scalar>* y,
  const Scalar alpha, 
  const Scalar beta) const 
{
  bool useA = alpha != Teuchos::ScalarTraits<Scalar>::one();
  bool useB = beta != Teuchos::ScalarTraits<Scalar>::zero();

  RCP<const Thyra::VectorBase<Scalar> > cxp = rcp(&x, false);
  RCP<Thyra::VectorBase<Scalar> > xp 
    = rcp_const_cast<Thyra::VectorBase<Scalar> >(cxp);
  Vector<Scalar> in = xp;
  Vector<Scalar> opOut;

  if (!useB)
  {
    /* wrap y in a temporary handle */
    opOut = rcp(y, false);
    applyOp(M_trans, in, opOut);
    if (useA)
    {
      opOut.scale(alpha);
    }
  }
  else
  {
    opOut = VectorSpace<Scalar>(this->range()).createMember();
    applyOp(M_trans, in, opOut);
    Vector<Scalar> Y = rcp(y, false);
    if (useA)
    {
      Y.update(alpha, opOut, beta);
    }
  }
}
    
 

/* ---- Simplified linear op with spaces ------- */

template <class Scalar> inline
SimplifiedLinearOpWithSpaces<Scalar>
::SimplifiedLinearOpWithSpaces(const VectorSpace<Scalar>& domain,
  const VectorSpace<Scalar>& range)
  : range_(range), domain_(domain) {}


template <class Scalar> inline
RCP< const VectorSpaceBase<Scalar> > 
SimplifiedLinearOpWithSpaces<Scalar>::range() const 
{
  return range_.ptr();
}


template <class Scalar> inline
RCP< const VectorSpaceBase<Scalar> > 
SimplifiedLinearOpWithSpaces<Scalar>::domain() const 
{
  return domain_.ptr();
}







}

#endif
