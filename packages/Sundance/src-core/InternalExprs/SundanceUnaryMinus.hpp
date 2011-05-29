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

#ifndef SUNDANCE_UNARYMINUS_H
#define SUNDANCE_UNARYMINUS_H

#include "SundanceDefs.hpp"
#include "SundanceUnaryExpr.hpp"
#include "SundanceUnaryMinusEvaluator.hpp"


namespace Sundance
{
  using namespace Sundance;
  using namespace Teuchos;
  

      /** */
      class UnaryMinus : public UnaryExpr,
                         GenericEvaluatorFactory<UnaryMinus, UnaryMinusEvaluator>
        {
        public:
          /** construct with the argument */
          UnaryMinus(const RCP<ScalarExpr>& arg);

          /** virtual dtor */
          virtual ~UnaryMinus() {;}


          /** */
          virtual std::ostream& toText(std::ostream& os, bool paren) const ;

          /** */
          virtual XMLObject toXML() const ;

          /** */
          virtual bool isLinear() const {return true;}

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

          /** */
          virtual Set<MultiSet<int> > internalFindQ_W(int order, 
                                                   const EvalContext& context) const ;

          /** */
          virtual RCP<ExprBase> getRcp() {return rcp(this);}

        };
    }

#endif
