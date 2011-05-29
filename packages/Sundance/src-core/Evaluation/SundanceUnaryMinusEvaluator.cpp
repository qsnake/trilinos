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

#include "SundanceUnaryMinusEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceUnaryMinus.hpp"

#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;





UnaryMinusEvaluator
::UnaryMinusEvaluator(const UnaryMinus* expr,
                      const EvalContext& context)
  : UnaryEvaluator<UnaryMinus>(expr, context)
{
  try
    {
      int vecResultIndex = 0;
      int constResultIndex = 0;
      
      for (int i=0; i<this->sparsity()->numDerivs(); i++)
        {
          /* Determine the index into which the result will be written */
          bool resultIsConstant = this->sparsity()->state(i)==ConstantDeriv; 
          
          if (!resultIsConstant)
            {
              addVectorIndex(i, vecResultIndex);
              vecResultIndex++;
            }
          else
            {
              addConstantIndex(i, constResultIndex);
              constResultIndex++;
            }
        }
    }
  catch(std::exception& e)
    {
      TEST_FOR_EXCEPTION(true, RuntimeError, 
                         "exception detected in UnaryMinusEvaluator: expr="
                         << expr->toString() << std::endl
                         << "exception=" << e.what());
    }
}

void UnaryMinusEvaluator
::internalEval(const EvalManager& mgr,
               Array<double>& constantResults,
               Array<RCP<EvalVector> >& vectorResults) const
{
  Tabs tab;
  SUNDANCE_MSG1(mgr.verb(),
    tab << "UnaryMinusEvaluator::eval() expr=" << expr()->toString());


  /* evaluate the argument */
  evalOperand(mgr, constantResults, vectorResults);


  if (mgr.verb() > 2)
    {
      Out::os() << tab << "UnaryMinus operand results" << std::endl;
      argSparsitySuperset()->print(Out::os(), vectorResults,
                           constantResults);
    }

  for (int i=0; i<constantResults.size(); i++)
    {
      constantResults[i] *= -1;
    }

  for (int i=0; i<vectorResults.size(); i++)
    {
      vectorResults[i]->multiply_S(-1.0);
    }

  
  if (mgr.verb() > 1)
    {
      Out::os() << tab << "UnaryMinus results" << std::endl;
      sparsity()->print(Out::os(), vectorResults,
                           constantResults);
    }

  
}


