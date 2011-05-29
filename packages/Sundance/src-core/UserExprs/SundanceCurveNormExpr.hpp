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

#ifndef SUNDANCE_CURVENORMEXPR_H
#define SUNDANCE_CURVENORMEXPR_H

#include "SundanceEvaluatorFactory.hpp"
#include "SundanceCurveNormEvaluator.hpp"

namespace Sundance
{
  using namespace Sundance;

  /** 
   * CurveNormExpr is an
   * expression that returns the Cartesian coordinate component
   * of the normal vector. Which dimensional component will be returned
   * can be specified in the constructor argument. <br>
   * The normal vector of the curve is always pointing outwards. This
   * outwards direction is in terms of the parameterized curve equation
   * from positive values to negative values.
   */
  class CurveNormExpr : public EvaluatableExpr,
                        public GenericEvaluatorFactory<CurveNormExpr, CurveNormEvaluator>
    {
    public:
      /** */
	  CurveNormExpr(int dir, const std::string& name="");

      /** */
      virtual ~CurveNormExpr() {;}

      /** */
      virtual XMLObject toXML() const ;

      /** */
      int dir() const {return dir_;}

      /** */
      const std::string& name() const {return name_;}

      /** Write a simple text description suitable 
       * for output to a terminal */
      virtual std::ostream& toText(std::ostream& os, bool paren) const 
        {os<<name(); return os;}
      
      /** */
      virtual Set<MultipleDeriv> 
      internalFindW(int order, const EvalContext& context) const ;
      
      /* we use the default implementation from the supper class */
      //virtual Set<MultipleDeriv>
      //internalFindC(int order, const EvalContext& context) const ;
      
      /*  we use the default implementation from the supper class */
      //virtual Set<MultipleDeriv>
      //internalFindV(int order, const EvalContext& context) const ;
      
      /** */
      virtual RCP<ExprBase> getRcp() {return rcp(this);}

      /** Ordering operator for use in transforming exprs to standard form */
      virtual bool lessThan(const ScalarExpr* other) const ;

      static std::string coordName(int dir, const std::string& name);

    private:
      int dir_;
      std::string name_;
    };
}

#endif
