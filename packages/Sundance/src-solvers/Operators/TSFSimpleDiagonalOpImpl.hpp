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

#ifndef TSF_SIMPLE_DIAGONAL_OP_IMPL_HPP
#define TSF_SIMPLE_DIAGONAL_OP_IMPL_HPP



#include "TSFSimpleDiagonalOpDecl.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorImpl.hpp"
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
SimpleDiagonalOp<Scalar>::SimpleDiagonalOp(
  const Vector<Scalar>& diag)
  : SimplifiedLinearOpWithSpaces<Scalar>(
    diag.space(), diag.space()
    ), diag_(diag)
{}
  
/* */
template <class Scalar> inline
void SimpleDiagonalOp<Scalar>::applyOp(const Thyra::EOpTransp M_trans,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab(0);
  SUNDANCE_MSG2(this->verb(), tab << "SimpleDiagonalOp::applyOp()");

  Vector<Scalar> tmp = in.dotStar(diag_);
  out.acceptCopyOf(tmp);

  SUNDANCE_MSG2(this->verb(), tab << "done SimpleDiagonalOp::applyOp()");
}
  
/* */
template <class Scalar> inline
std::string SimpleDiagonalOp<Scalar>::description() const 
{
  return "DiagonalOp[diag=" + diag_.description() + "]";
}


/* */
template <class Scalar> inline
void SimpleDiagonalOp<Scalar>::print(std::ostream& os) const 
{
  Tabs tab(0);
  os << tab << "DiagonalOp[" << std::endl;
  Tabs tab1;
  os << tab1 << "diag = " << diag_ << std::endl;
  os << tab << "]" << std::endl;
}



template <class Scalar> inline
LinearOperator<Scalar> diagonalOperator(
  const Vector<Scalar>& diag)
{
  RCP<LinearOpBase<Scalar> > A 
    = rcp(new SimpleDiagonalOp<Scalar>(diag));

  return A;
}



}

#endif
