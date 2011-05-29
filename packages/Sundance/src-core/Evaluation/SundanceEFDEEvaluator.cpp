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

#include "SundanceEFDEEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceSet.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;


EFDEEvaluator::EFDEEvaluator(
  const ExplicitFunctionalDerivativeElement* expr, 
  const EvalContext& context
  )
  : UnaryEvaluator<ExplicitFunctionalDerivativeElement>(expr, context),
    constValIndexToArgIndexMap_(),
    varValIndexToArgIndexMap_()
{

  Tabs tabs;
  SUNDANCE_VERB_LOW(tabs << "initializing EFDE evaluator for " 
                    << expr->toString());
  SUNDANCE_VERB_MEDIUM(tabs << "return sparsity " << std::endl << *(this->sparsity)());

  /* 
   * This evaluator requires no calculations. All that is done is to
   * map derivatives (md, fd) in the argument's result arrays to 
   * derivatives (md) in this expression's result arrays. 
   */
  

  int vecResultIndex = 0;
  int constResultIndex = 0;

  const Deriv& fd = expr->fd();

  for (int i=0; i<this->sparsity()->numDerivs(); i++)
    {
      const MultipleDeriv& d = this->sparsity()->deriv(i);
      const DerivState& myState = this->sparsity()->state(i);

      if (myState==ConstantDeriv)
      {
        Tabs tab2;
        SUNDANCE_VERB_HIGH(tab2 
          << "deriv is constant, will be stored at index "
          << constResultIndex << " in the const result array");
        addConstantIndex(i, constResultIndex++);
      }
      else
      {
        Tabs tab2;
        SUNDANCE_VERB_HIGH(tab2 
          << "deriv is variable, will be stored at index "
          << vecResultIndex << " in the var result array");
        addVectorIndex(i, vecResultIndex++);
      }
      
      MultipleDeriv dArg = d;
      dArg.put(fd);

      int argIndex = argSparsitySuperset()->getIndex(dArg);

      
      TEST_FOR_EXCEPTION(argIndex==-1, RuntimeError,
        "Derivative " << dArg << " expected in argument but not found");

      
      const DerivState& argState = argSparsitySuperset()->state(argIndex);
      TEST_FOR_EXCEPTION(argState != myState, InternalError, 
        "mismatched states");

      if (argState==ConstantDeriv)
      {
        int constArgIndex = argEval()->constantIndexMap().get(argIndex);
        constValIndexToArgIndexMap_.append(constArgIndex);
      }
      else
      {
        int vectorArgIndex = argEval()->vectorIndexMap().get(argIndex);
        varValIndexToArgIndexMap_.append(vectorArgIndex);
      }
    }
  
  SUNDANCE_VERB_HIGH(tabs 
    << " constant index map " 
    << constValIndexToArgIndexMap_ << std::endl 
    << " vector index map " 
    << varValIndexToArgIndexMap_
    );
}



void EFDEEvaluator::internalEval(const EvalManager& mgr,
  Array<double>& constantResults,
  Array<RCP<EvalVector> >& vectorResults) const 
{
  TimeMonitor timer(efdeEvalTimer());
  Tabs tabs;

  SUNDANCE_MSG1(mgr.verb(), tabs << "EFDEEvaluator::eval() expr=" 
    << expr()->toString());

  SUNDANCE_MSG2(mgr.verb(), tabs << "sparsity = " << std::endl 
    << *(this->sparsity)());

  constantResults.resize(constValIndexToArgIndexMap_.size());
  vectorResults.resize(varValIndexToArgIndexMap_.size());

  /* evaluate the argument */
  Array<RCP<EvalVector> > argVectorResults;
  Array<double> argConstantResults;

  evalOperand(mgr, argConstantResults, argVectorResults);


  if (mgr.verb() > 2)
    {
      Tabs tab1;
      Out::os() << tab1 << "EFDE operand results" << std::endl;
      argSparsitySuperset()->print(Out::os(), argVectorResults,
                                   argConstantResults);
    }


  for (int i=0; i<constantResults.size(); i++)
  {
    constantResults[i] = argConstantResults[constValIndexToArgIndexMap_[i]];
  }

  
  for (int i=0; i<vectorResults.size(); i++)
  {
    vectorResults[i] = mgr.popVector();
    const RCP<EvalVector>& v = argVectorResults[varValIndexToArgIndexMap_[i]];
    vectorResults[i]->setTo_V(v.get());
  }

  
  

  if (mgr.verb() > 2)
  {
    Tabs tab1;
    Out::os() << tab1 << "results " << std::endl;
    this->sparsity()->print(Out::os(), vectorResults,
      constantResults);
  }
}

