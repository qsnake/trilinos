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

#include "SundanceEvaluatorFactory.hpp"
#include "SundanceInstructionCachingEvaluator.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceCurveNormExpr.hpp"
#include "SundanceCurveNormEvaluator.hpp"
#include "SundanceSpatiallyConstantExpr.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceDiscreteFuncElement.hpp"

using namespace Sundance;
using namespace Teuchos;


EvaluatorFactory::EvaluatorFactory()
{;}

Evaluator* EvaluatorFactory::commonCreate(const EvaluatableExpr* expr,
                                          const EvalContext& context,
                                          int topLevelDiffOrder) const
{
  const CoordExpr* c = dynamic_cast<const CoordExpr*>(expr);
  if (c != 0)
    {
      return new CoordExprEvaluator(c, context, topLevelDiffOrder);
    }
  
  const SpatiallyConstantExpr* sc 
    = dynamic_cast<const SpatiallyConstantExpr*>(expr);
  if (sc != 0)
    {
      return new ConstantEvaluator(sc, context, topLevelDiffOrder);
    }

  const SymbolicFuncElement* u
    = dynamic_cast<const SymbolicFuncElement*>(expr);
  if (u != 0)
    {
      return new SymbolicFuncElementEvaluator(u, context, topLevelDiffOrder);
    }
  
  const DiscreteFuncElement* df
    = dynamic_cast<const DiscreteFuncElement*>(expr);
  if (df != 0)
    {
      return new DiscreteFuncElementEvaluator(df, context, topLevelDiffOrder);
    }

  const CurveNormExpr* cne
    = dynamic_cast<const DiscreteFuncElement*>(expr);
  if (cne != 0)
    {
      return new CurveNormEvaluator(cne, context, topLevelDiffOrder);
    }

  TEST_FOR_EXCEPTION(true, InternalError,
                     "EvaluatorFactory::commonCreate() could not create an "
                     "evaluator for " << expr->toString());

  return 0;
}


RCP<EvaluatorFactory>&  EvaluatorFactory::defaultEvaluator()
{
  static RCP<EvaluatorFactory> rtn 
    = rcp(new InstructionCachingEvaluatorFactory());
  return rtn;
}
