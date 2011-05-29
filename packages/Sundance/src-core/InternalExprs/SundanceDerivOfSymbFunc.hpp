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

#ifndef SUNDANCE_DERIVOFSYMBFUNC_H
#define SUNDANCE_DERIVOFSYMBFUNC_H

#include "SundanceDefs.hpp"
#include "SundanceFunctionIdentifier.hpp"
#include "SundanceDiffOp.hpp"
#include "SundanceDerivOfSymbFuncEvaluator.hpp"


namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;




/**
 * Specialization of DiffOp to the case where the argument is a 
 * symbolic function, allowing optimized evaluation.
 */
class DerivOfSymbFunc : public DiffOp
{
public:
  /** ctor */
  DerivOfSymbFunc(const MultiIndex& op, 
    const RCP<ScalarExpr>& arg);

  /** virtual destructor */
  virtual ~DerivOfSymbFunc() {;}

  /** */
  const FunctionIdentifier& argFid() const {return argFid_;}

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}

  /** */
  virtual Evaluator* createEvaluator(const EvaluatableExpr* expr,
    const EvalContext& context) const ;

  /** */
  Deriv representMeAsFunctionalDeriv() const ;

  /** Ordering operator for use in transforming exprs to standard form */
  virtual bool lessThan(const ScalarExpr* other) const ;


private:
  FunctionIdentifier argFid_;
};
}

#endif
