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

#ifndef SUNDANCE_SCALAREXPR_H
#define SUNDANCE_SCALAREXPR_H


#include "SundanceDefs.hpp"
#include "SundanceExceptions.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_RefCountPtrDecl.hpp"
#include "SundanceExprBase.hpp"
#include "SundanceExpr.hpp"

namespace Sundance
{

using namespace Sundance;
using namespace Teuchos;





/** */
class ScalarExpr : virtual public ExprBase
{
public:
  /** empty ctor */
  ScalarExpr();

  /** virtual destructor */
  virtual ~ScalarExpr() {;}


  /** Indicate whether this expression is constant in space */
  virtual bool isConstant() const {return false;}


  /** Indicate whether this expression is an immutable constant */
  virtual bool isImmutable() const {return false;}

  /** Indicate whether this expression is a "hungry"
   * differential operator that is awaiting an argument. */
  virtual bool isHungryDiffOp() const {return false;}

  /** Indicate whether the expression is independent of the given 
   * functions */
  virtual bool isIndependentOf(const Expr& u) const {return true;}

  /** 
   * Indicate whether the expression is nonlinear 
   * with respect to test functions */
  virtual bool isLinearInTests() const {return false;}

  /** 
   * Indicate whether every term in the expression contains test functions */
  virtual bool everyTermHasTestFunctions() const {return hasTestFunctions();}

  /** 
   * Indicate whether the expression contains test functions */
  virtual bool hasTestFunctions() const {return false;}
  /** 
   * Indicate whether the expression contains unknown functions */
  virtual bool hasUnkFunctions() const {return false;}

  /** Indicate whether the expression is linear in the given 
   * functions */
  virtual bool isLinearForm(const Expr& u) const {return false;}

  /** Indicate whether the expression is quadratic in the given 
   * functions */
  virtual bool isQuadraticForm(const Expr& u) const {return false;}

  /** Find all the unknown functions in this expression. */
  virtual void getUnknowns(Set<int>& unkID, Array<Expr>& unks) const {;}

  /** Find all the test functions in this expression. */
  virtual void getTests(Set<int>& varID, Array<Expr>& vars) const {;}

  /** Ordering operator for use in transforming exprs 
   * to standard form */
  virtual bool lessThan(const ScalarExpr* other) const = 0 ;


protected:
private:
};
}

#endif
