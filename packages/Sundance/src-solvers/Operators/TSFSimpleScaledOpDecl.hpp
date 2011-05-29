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

#ifndef TSF_SIMPLE_SCALED_OP_DECL_HPP
#define TSF_SIMPLE_SCALED_OP_DECL_HPP



#include "SundanceDefs.hpp"
#include "TSFSimplifiedLinearOpBaseDecl.hpp"
#include "TSFLinearOperatorDecl.hpp"


namespace TSFExtended
{
using namespace Teuchos;
using namespace Sundance;



/**
 * Represent a scaled operator
 */
template <class Scalar>
class SimpleScaledOp : public SimplifiedLinearOpWithSpaces<Scalar>,
                       public Printable
{
public:
  /** */
  SimpleScaledOp(const Scalar& alpha, const LinearOperator<Scalar>& A);
  
  /** */
  void applyOp(const Thyra::EOpTransp M_trans,
    const Vector<Scalar>& in,
    Vector<Scalar> out) const;

  
  /** */
  std::string description() const ;

  /** */
  const Scalar& alpha() const {return alpha_;}

  /** */
  LinearOperator<Scalar> op() const {return A_;}

  /** */
  void print(std::ostream& os) const ;

private:
  Scalar alpha_;
  LinearOperator<Scalar> A_;
};


/** \relates SimpleScaledOp */
template <class Scalar>
LinearOperator<Scalar> scaledOperator(
  const Scalar& scale,
  const LinearOperator<Scalar>& op);


template <class Scalar> 
LinearOperator<Scalar> operator*(const Scalar& a, const LinearOperator<Scalar>& A);

}

#endif
