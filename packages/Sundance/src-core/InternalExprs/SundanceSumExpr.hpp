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

#ifndef SUNDANCE_SUMEXPR_H
#define SUNDANCE_SUMEXPR_H

#include "SundanceBinaryExpr.hpp"
#include "SundanceSumEvaluator.hpp"

namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;



/**
 * SumExpr is the internal representation of an addition or subtraction
 * node in the expression tree.
 */
class SumExpr : public BinaryExpr,
                public GenericEvaluatorFactory<SumExpr, SumEvaluator>
{
public:
  /** */
  SumExpr(const RCP<ScalarExpr>& a, 
    const RCP<ScalarExpr>& b, int sign);

  /** virtual dtor */
  virtual ~SumExpr() {;}

  /** */
  virtual bool isHungryDiffOp() const ;

  /** */
  virtual bool isLinear() const {return true;}

          

  /** 
   * Indicate whether the expression is linear 
   * with respect to test functions */
  virtual bool isLinearInTests() const ;

  /** 
   * Indicate whether every term in the expression contains test functions */
  virtual bool everyTermHasTestFunctions() const ;
          

  /** Indicate whether the expression is linear in the given 
   * functions */
  virtual bool isLinearForm(const Expr& u) const ;

  /** Indicate whether the expression is a 
   * quadratic form in the given functions */
  virtual bool isQuadraticForm(const Expr& u) const ;

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}
          
  /** */
  const Map<Expr, int>& getSumTree() const {return sumTree_;}

protected:
  /** */
  virtual bool parenthesizeSelf() const {return true;}
  /** */
  virtual bool parenthesizeOperands() const {return false;}
  /** */
  virtual const std::string& xmlTag() const ;
  /** */
  virtual const std::string& opChar() const ;

private:
  Map<Expr, int> sumTree_;


};
}
#endif
