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

#ifndef SUNDANCE_DIFFOP_H
#define SUNDANCE_DIFFOP_H

#include "SundanceDefs.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnaryExpr.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceMap.hpp"
#include "SundanceSet.hpp"
#include "SundanceMultipleDeriv.hpp"
#include "SundanceDiffOpEvaluator.hpp"


namespace Sundance
{
class Expr;

using namespace Sundance;
using namespace Teuchos;




/**
 *
 */
class DiffOp : public UnaryExpr
{
public:
  /** ctor */
  DiffOp(const MultiIndex& op, const RCP<ScalarExpr>& arg);

  /** virtual destructor */
  virtual ~DiffOp() {;}


  /** 
   * Indicate whether the expression is linear 
   * with respect to test functions 
   */
  virtual bool isLinearInTests() const 
    {return evaluatableArg()->isLinearInTests();}
      
      
  /** Indicate whether the expression is linear in the given 
   * functions */
  virtual bool isLinearForm(const Expr& u) const 
    {
      return evaluatableArg()->isLinearForm(u);
    }
      
  /** Indicate whether the expression is at most 
   * quadratic in the given functions */
  virtual bool isQuadraticForm(const Expr& u) const
    {
      return evaluatableArg()->isQuadraticForm(u);
    }

  /**
   * Find the maximum differentiation order acting on discrete
   * functions in this expression. 
   */
  virtual int maxDiffOrderOnDiscreteFunctions() const 
    {
      int rtn = evaluatableArg()->maxDiffOrderOnDiscreteFunctions();
      if (evaluatableArg()->hasDiscreteFunctions()) 
      {
        rtn += mi_.order();
      }
      return rtn;
    }


  /** Write a simple text description suitable
   * for output to a terminal */
  virtual std::ostream& toText(std::ostream& os, bool paren) const ;

  /** Write in XML */
  virtual XMLObject toXML() const ;


  /** */
  virtual Set<MultipleDeriv> 
  internalFindW(int order, const EvalContext& context) const ;


  /** */
  virtual Set<MultipleDeriv> 
  internalFindV(int order, const EvalContext& context) const ;


  /** */
  virtual Set<MultipleDeriv> 
  internalFindC(int order, const EvalContext& context) const ;


  /** */
  virtual RCP<Array<Set<MultipleDeriv> > > 
  internalDetermineR(const EvalContext& context,
    const Array<Set<MultipleDeriv> >& RInput) const ;

  /** */
  void requestMultiIndexAtEvalPoint(const MultiIndex& mi,
    const MultipleDeriv& u,
    const EvalContext& context) const ;
      
      
      
  /** */
  const Deriv& myCoordDeriv() const {return myCoordDeriv_;}

  /** */
  const MultiIndex& mi() const {return mi_;}

  /** Get the functions that are required in the evaluation
   * of the multiple deriv d */
  const Sundance::Set<Deriv>& requiredFunctions(const MultipleDeriv& d) const 
    {return requiredFunctions_[d];}

  /** */
  bool requiresFunctionsToEval(const MultipleDeriv& d) const 
    {return requiredFunctions_.containsKey(d);}

    
      
     

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}


  /** */
  virtual Evaluator* createEvaluator(const EvaluatableExpr* expr,
    const EvalContext& context) const ;

     

  /** */
  virtual void registerSpatialDerivs(const EvalContext& context, 
    const Set<MultiIndex>& miSet) const ;

  /** Ordering operator for use in transforming exprs to standard form */
  virtual bool lessThan(const ScalarExpr* other) const ;

private:

      
      

  MultiIndex mi_;


  Deriv myCoordDeriv_;

  mutable Map<MultipleDeriv, Sundance::Set<Deriv>, 
              increasingOrder<MultipleDeriv> > requiredFunctions_;

  mutable bool ignoreFuncTerms_;
};
}


#endif
