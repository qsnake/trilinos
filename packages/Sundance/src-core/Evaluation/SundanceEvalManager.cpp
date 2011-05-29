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


#include "SundanceEvalManager.hpp"
#include "SundanceOut.hpp"
#include "SundanceAbstractEvalMediator.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceCellDiameterExpr.hpp"
#include "SundanceCurveNormExpr.hpp"
#include "SundanceCellVectorExpr.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceMultiIndex.hpp"



using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

EvalManager::EvalManager()
  : verb_(0),
    region_(),
    mediator_()
{}

void EvalManager::setVerbosity(int verb)
{
  verb_ = verb;
//  if (mediator_.get()) mediator_->setVerbosity(verb);
}

void EvalManager::evalCoordExpr(const CoordExpr* expr,
                                RCP<EvalVector>& result) const 
{

  TimeMonitor timer(coordEvalTimer());
  TEST_FOR_EXCEPTION(mediator() == 0, InternalError,
                     "uninitialized mediator in "
                     "EvalManager::evalCoordExpr");
  mediator()->evalCoordExpr(expr, result);
}


void EvalManager::evalCellDiameterExpr(const CellDiameterExpr* expr,
                                RCP<EvalVector>& result) const 
{
  TEST_FOR_EXCEPTION(mediator() == 0, InternalError,
                     "uninitialized mediator in "
                     "EvalManager::evalCellDiameterExpr");
  mediator()->evalCellDiameterExpr(expr, result);
}

void EvalManager::evalCurveNormExpr(const CurveNormExpr* expr,
                                RCP<EvalVector>& result) const
{
  TEST_FOR_EXCEPTION(mediator() == 0, InternalError,
                     "uninitialized mediator in "
                     "EvalManager::evalCurveNormExpr");
  mediator()->evalCurveNormExpr(expr, result);
}

void EvalManager::evalCellVectorExpr(const CellVectorExpr* expr,
                                RCP<EvalVector>& result) const 
{
  TEST_FOR_EXCEPTION(mediator() == 0, InternalError,
                     "uninitialized mediator in "
                     "EvalManager::evalCellVectorExpr");
  mediator()->evalCellVectorExpr(expr, result);
}


void EvalManager::evalDiscreteFuncElement(const DiscreteFuncElement* expr,
                                          const Array<MultiIndex>& mi,
                                          Array<RCP<EvalVector> >& result) const 
{
  TimeMonitor timer(discFuncEvalTimer());

  TEST_FOR_EXCEPTION(mediator() == 0, InternalError,
                     "uninitialized mediator in "
                     "EvalManager::evalDiscreteFuncElement");

  mediator()->evalDiscreteFuncElement(expr, mi, result);
}


RCP<EvalVector> EvalManager::popVector() const
{
  return stack().popVector();
}

TempStack& EvalManager::stack()
{
  static TempStack rtn;
  return rtn;
}
