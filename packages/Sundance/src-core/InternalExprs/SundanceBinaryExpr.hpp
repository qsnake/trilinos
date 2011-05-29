/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
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
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#ifndef SUNDANCE_BINARYEXPR_H
#define SUNDANCE_BINARYEXPR_H

#include "SundanceDefs.hpp"
#include "SundanceExprWithChildren.hpp"
#include "SundanceExpr.hpp"

namespace Sundance
{
using namespace Teuchos;

/** 
 * BinaryExpr is a base class for binary expressions, e.g., sums
 * and products. It provides a number of helper methods.
 */
class BinaryExpr : public ExprWithChildren
{
public:
  /** construct with left and right operands */
  BinaryExpr(const RCP<ScalarExpr>& left,
    const RCP<ScalarExpr>& right, int sign);

  /** virtual dtor */
  virtual ~BinaryExpr() {;}

  /** */
  virtual std::ostream& toText(std::ostream& os, bool paren) const ;

  /** */
  virtual XMLObject toXML() const ;

  /** */
  Expr left() const {return child(0);}

  /** */
  Expr right() const {return child(1);}

  /** */
  int sign() const {return sign_;}

  /** Downcast the left expr to an evaluatable expr */
  const EvaluatableExpr* leftEvaluatable() const 
    {return evaluatableChild(0);}

  /** Downcast the right expr to an evaluatable expr */
  const EvaluatableExpr* rightEvaluatable() const 
    {return evaluatableChild(1);}

  /** Downcast the left expr to a scalar expr */
  const ScalarExpr* leftScalar() const {return scalarChild(0);}

  /** Downcast the right expr to a scalar expr */
  const ScalarExpr* rightScalar() const {return scalarChild(1);}

  /** Ordering operator for use in transforming exprs to standard form */
  virtual bool lessThan(const ScalarExpr* other) const ;

protected:

          

  /** */
  virtual bool parenthesizeSelf() const = 0 ;
  /** */
  virtual bool parenthesizeOperands() const = 0 ;
  /** */
  virtual const std::string& xmlTag() const = 0 ;
  /** */
  virtual const std::string& opChar() const = 0 ;



private:
          

  int sign_;
};
}

#endif
