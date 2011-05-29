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

#ifndef SUNDANCE_CONSTANTEXPR_H
#define SUNDANCE_CONSTANTEXPR_H

#include "SundanceSpatiallyConstantExpr.hpp"

namespace Sundance
{

/**
 * ConstantExpr contains an immutable constant, to be distinguished
 * from a parameter that is constant in space but can change
 * during the course of a simulation.
 */
class ConstantExpr : public SpatiallyConstantExpr
{
public:
  ConstantExpr(const double& value);
  virtual ~ConstantExpr() {;}

  /** */
  virtual std::ostream& toText(std::ostream& os, bool paren) const ;

  /** */
  virtual XMLObject toXML() const ;

  /** */
  virtual bool isImmutable() const {return true;}

          
  /** */
  virtual void setValue(const double& value) {value_ = value;}
          
  /** */
  virtual const double& value() const {return value_;}

  /** Ordering operator for use in transforming exprs to standard form */
  virtual bool lessThan(const ScalarExpr* other) const ;

  /** */
  virtual Set<MultipleDeriv> 
  internalFindW(int order, const EvalContext& context) const ;

  /** Find spatially-constant functional derivatives */
  virtual Set<MultipleDeriv> 
  internalFindC(int order, const EvalContext& context) const ;

  /** Find spatially-variable functional derivatives */
  virtual Set<MultipleDeriv> 
  internalFindV(int order, const EvalContext& context) const ;
          


  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}


protected:
private:
  double value_;
};
}

#endif
