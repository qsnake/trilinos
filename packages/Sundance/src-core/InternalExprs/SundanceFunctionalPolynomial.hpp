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

#ifndef SUNDANCE_FUNCTIONALPOLYNOMIAL_H
#define SUNDANCE_FUNCTIONALPOLYNOMIAL_H

#include "SundanceDefs.hpp"
#include "SundanceEvaluatableExpr.hpp"


namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;




/**
 * Specialized class for representing polynomials in symbolic
 * functions and their derivatives. 
 */
class FunctionalPolynomial : public EvaluatableExpr
{
public:
  /** ctor */
  FunctionalPolynomial(const RCP<ScalarExpr>& expr);
  /** ctor */
  FunctionalPolynomial(const Map<int, RCP<ScalarExpr> >& funcs,
    const Map<int, Set<MultiIndex> >& funcMultiIndices,
    const Array<Map<MultipleDeriv, RCP<ScalarExpr> > > & coeffs);

  /** virtual destructor */
  virtual ~FunctionalPolynomial() {;}


  /** */
  virtual Set<MultipleDeriv> 
  internalFindW(int order, const EvalContext& context) const ;

      

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}

  /** */
  virtual Evaluator* createEvaluator(const EvaluatableExpr* expr,
    const EvalContext& context) const ;

  /** */
  RCP<FunctionalPolynomial> addPoly(const FunctionalPolynomial* other,
    int sign) const ;

  /** */
  RCP<FunctionalPolynomial> multiplyPoly(const FunctionalPolynomial* other) const ;

  /** */
  RCP<FunctionalPolynomial> multiplyScalar(const RCP<ScalarExpr>& alpha) const ;

  /** */
  RCP<FunctionalPolynomial> addFunction(const RCP<ScalarExpr>& u,
    int sign) const ;

  /** */
  static bool isConvertibleToPoly(const ScalarExpr* expr) ;

  /** */
  static RCP<FunctionalPolynomial> toPoly(const RCP<ScalarExpr>& expr);


  /** Write a simple text description suitable 
   * for output to a terminal */
  virtual std::ostream& toText(std::ostream& os, bool paren) const ;
      
  /** Write in XML */
  virtual XMLObject toXML() const ;

  /** Ordering operator for use in transforming exprs to standard form */
  virtual bool lessThan(const ScalarExpr* other) const ;

private:

  /** */
  Map<int, RCP<ScalarExpr> > funcs_;

  /** */
  Map<int, Set<MultiIndex> > funcMultiIndices_;

  /** */
  Array<Map<MultipleDeriv, RCP<ScalarExpr> > > coeffs_;

  /** */
  Array<Set<MultipleDeriv> > keys_;

  /** */
  Set<Deriv> findFuncsForSummation(const Set<MultipleDeriv>& prevSet,
    const MultipleDeriv& thisSet) const ;
      
  /**    
   * Given a term's key, find the term that will appear as
   * its successor in the evaluation recurrence.
   */
  MultipleDeriv successorTerm(const MultipleDeriv& md) const ;

  /** */
  void stepRecurrence(int level, const Map<MultipleDeriv, std::string>& sPrev,
    Map<MultipleDeriv, std::string>& sCurr) const ;

  /** */
  std::string evalString() const ;
};
}

#endif
