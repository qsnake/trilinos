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



#ifndef SUNDANCE_PRODUCTEXPR_H
#define SUNDANCE_PRODUCTEXPR_H

#include "SundanceBinaryExpr.hpp"


namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;





/** 
 * ProductExpr represents a product of two scalar-valued expressions
 */
class ProductExpr : public BinaryExpr
{
public:
  /** */
  ProductExpr(const RCP<ScalarExpr>& a, 
    const RCP<ScalarExpr>& b);

  /** virtual dtor */
  virtual ~ProductExpr() {;}

  /** Indicate whether this expression is a "hungry"
   * differential operator that is awaiting an argument. */
  virtual bool isHungryDiffOp() const ;

  /** */
  virtual Evaluator* createEvaluator(const EvaluatableExpr* expr,
    const EvalContext& context) const ;
      
  /** */
  virtual Set<MultiSet<int> > internalFindQ_W(int order, 
    const EvalContext& context) const ;
      
  /** */
  virtual Set<MultiSet<int> > internalFindQ_V(int order, 
    const EvalContext& context) const ;

      

  /** */
  virtual bool isProduct() const {return true;}

  /** 
   * Indicate whether the expression is linear 
   * with respect to test functions */
  virtual bool isLinearInTests() const 
    {
      bool leftIsLinear = leftScalar()->isLinearInTests() ;
      bool rightIsLinear = rightScalar()->isLinearInTests() ;

      bool leftHasTests = leftScalar()->hasTestFunctions() ;
      bool rightHasTests = rightScalar()->hasTestFunctions() ;
          
      return (leftIsLinear && !rightHasTests) 
        || (!leftHasTests && rightIsLinear);
    }
      
  /** Indicate whether the expression is linear in the given 
   * functions */
  virtual bool isLinearForm(const Expr& u) const 
    {
      bool L = leftScalar()->isLinearForm(u);
      bool R = rightScalar()->isLinearForm(u);
      return ( (L || R) && !(L==R) );
    }
      
  /** Indicate whether the expression is a quadratic form in the given 
   * functions */
  virtual bool isQuadraticForm(const Expr& u) const
    {
      bool LQ = leftScalar()->isQuadraticForm(u);
      bool RQ = rightScalar()->isQuadraticForm(u);
      bool LL = leftScalar()->isLinearForm(u);
      bool RL = leftScalar()->isLinearForm(u);
      bool LI = leftScalar()->isIndependentOf(u);
      bool RI = rightScalar()->isIndependentOf(u);

      if (LI && RQ) return true;
      if (RI && LQ) return true;
      if (LL && RL) return true;
      return false;
    }

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}
protected:
  /** */
  virtual bool parenthesizeSelf() const {return true;}
  /** */
  virtual bool parenthesizeOperands() const {return true;}
  /** */
  virtual const std::string& xmlTag() const ;
  /** */
  virtual const std::string& opChar() const ;

private:

};
}

#endif
