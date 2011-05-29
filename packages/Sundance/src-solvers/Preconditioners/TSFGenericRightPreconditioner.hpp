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

#ifndef TSFGENERICRIGHTPRECONDITIONER_HPP
#define TSFGENERICRIGHTPRECONDITIONER_HPP

#include "SundanceDefs.hpp"
#include "TSFVectorDecl.hpp"
#include "TSFLinearOperatorDecl.hpp"
#include "TSFPreconditionerBase.hpp"


namespace TSFExtended
{
using namespace Teuchos;

/**
 * A one-size-fits-most right preconditioner that can be constructed by
 * accepting an operator for the right op of the preconditioner. 
 */
template <class Scalar>
class GenericRightPreconditioner : public PreconditionerBase<Scalar>
{
public:
  /** construct with an operator for the right preconditioner */
  GenericRightPreconditioner(const LinearOperator<Scalar>& right) 
    : PreconditionerBase<Scalar>(), right_(right) {;}

  /** virtual dtor */
  virtual ~GenericRightPreconditioner(){;}

    
  /** Return the right operator */
  virtual LinearOperator<Scalar> right() const {return right_;}

  /** A call to left() results in an error for a right precond. */
  virtual LinearOperator<Scalar> left() const
    {
      TEST_FOR_EXCEPTION(true, std::logic_error, "left() called for a "
        "preconditioner known to be a right precond");
      return LinearOperator<Scalar>();
    }

  /** return true because 
   * this preconditioner has a nontrivial right component. */
  virtual bool hasRight() const {return true;}

  /** return false, because this preconditioner has
   * no nontrivial left component */
  virtual bool hasLeft() const {return false;}

  /* Handleable boilerplate */
  GET_RCP(PreconditionerBase<Scalar>);

private:
  LinearOperator<Scalar> right_;
};



}

#endif
