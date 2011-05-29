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

#ifndef SUNDANCE_SYMBOLICFUNCEVALUATOR_H
#define SUNDANCE_SYMBOLICFUNCEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceSubtypeEvaluator.hpp"
#include "Teuchos_TimeMonitor.hpp"


namespace Sundance 
{
class Parameter;

class DiscreteFuncElementEvaluator; 
class SymbolicFuncElement;
class ConstantEvaluator; 

/** 
 *
 */
class SymbolicFuncElementEvaluator 
  : public SubtypeEvaluator<SymbolicFuncElement>
{
public:
  /** */
  SymbolicFuncElementEvaluator(const SymbolicFuncElement* expr, 
    const EvalContext& context);

  /** */
  virtual ~SymbolicFuncElementEvaluator(){;}

  /** */
  virtual void internalEval(const EvalManager& mgr,
    Array<double>& constantResults,
    Array<RCP<EvalVector> >& vectorResults) const ;

  /** */
  TEUCHOS_TIMER(symbolicFuncEvalTimer, "symbolic function evaluation");

  /** */
  const DiscreteFuncElementEvaluator* dfEval() const {return dfEval_;}
  /** */
  const ConstantEvaluator* pEval() const {return pEval_;}

private:
  Array<MultiIndex> mi_;
  Array<int> spatialDerivPtrs_;
  Array<int> onePtrs_;
  Array<int> paramValuePtrs_;
  const DiscreteFuncElement* df_;
  const Parameter* p_;
  const DiscreteFuncElementEvaluator* dfEval_;
  const ConstantEvaluator* pEval_;
  Array<string> stringReps_;
};

    

}

#endif
