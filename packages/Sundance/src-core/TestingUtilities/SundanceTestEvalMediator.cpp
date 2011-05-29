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


#include "SundanceTestEvalMediator.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceCellDiameterExpr.hpp"
#include "SundanceCellVectorExpr.hpp"
#include "SundanceDiscreteFuncElement.hpp"

#include "SundanceEvalManager.hpp"
#include "SundanceOut.hpp"
#include "SundanceExpr.hpp"
#include "SundanceTabs.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceEvalManager.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace SundanceTesting;

using namespace Teuchos;
using namespace std;



TestEvalMediator::TestEvalMediator(const Expr& fields)
  : AbstractEvalMediator(),
    x_(),
    funcIdToFieldNumberMap_(),
    fields_(fields.totalSize()),
    fieldNames_(fields.totalSize())
{
  EvalManager::stack().setVecSize(1);

  Expr f = fields.flatten();
  for (int i=0; i<f.size(); i++)
    {
      const DiscreteFuncElement* u0 
        = dynamic_cast<const DiscreteFuncElement*>(f[i].ptr().get());
      TEST_FOR_EXCEPTION(u0 == 0, InternalError,
                         "TestEvalMediator ctor: field argument "
                         << f[i] << " is not a discrete function");
      funcIdToFieldNumberMap_.put(u0->fid().dofID(), i);

      RCP<const DiscreteFuncDataStub> data = u0->commonData();
      const TestDiscreteFuncData* tdfd  
        = dynamic_cast<const TestDiscreteFuncData*>(data.get());

      TEST_FOR_EXCEPTION(tdfd==0, InternalError,
                         "df " << f[i] << " is not a TestDiscreteFunction");
      TEST_FOR_EXCEPTION(tdfd==0, InternalError,
                         "TestEvalMediator ctor: field argument "
                         << f[i] << " is not a TestDiscreteFunction");
      fields_[i] = tdfd->field();
      fieldNames_[i] = f[i].toString();
    }
}



void TestEvalMediator::evalCoordExpr(const CoordExpr* expr,
                                     RCP<EvalVector>& vec) const
{
  Tabs tabs;
  SUNDANCE_MSG1(verb(), tabs << "evaluating coord expr " << expr->toXML().toString());
  
  vec->setString(expr->name());

  int direction = expr->dir();
  
  double * const xx = vec->start();

  xx[0] = x_[direction];

  SUNDANCE_MSG2(verb(), tabs << "results: " << *vec);
}

void TestEvalMediator::evalCellDiameterExpr(const CellDiameterExpr* expr,
                                     RCP<EvalVector>& vec) const
{
  Tabs tabs;

  SUNDANCE_MSG1(verb(),
               tabs << "evaluating cell diameter expr " << expr->toXML().toString());
  
  vec->setString(expr->name());

  double * const xx = vec->start();

  xx[0] = 1.0;

  SUNDANCE_MSG2(verb(), tabs << "results: " << *vec);
}


void TestEvalMediator::evalCellVectorExpr(const CellVectorExpr* expr,
                                     RCP<EvalVector>& vec) const
{
  Tabs tabs;

  SUNDANCE_MSG1(verb(),
               tabs << "evaluating cell vector expr " << expr->toXML().toString());
  
  vec->setString(expr->name());

  int dim = expr->dimension();
  double * const xx = vec->start();

  if (expr->isNormal())
    {
      int c = expr->componentIndex();
      if (dim==1)
	{
	  xx[0] = 1.0;
	}
      else if (dim==2)
	{
	  if (c==0) xx[0] = 0.5;
	  else xx[0] = ::sqrt(3.0)/2.0;
	}
      else 
	{
	  if (c==0) xx[0] = 0.5;
	  else if (c==1) xx[0] = ::sqrt(3.0)/2.0 * 0.5;
	  else xx[0] = ::sqrt(3.0)/2.0 * ::sqrt(3.0)/2.0;
	}
    }
  TEST_FOR_EXCEPT(expr->isTangent());
  SUNDANCE_MSG2(verb(), tabs << "results: " << *vec);
}



void TestEvalMediator
::evalDiscreteFuncElement(const DiscreteFuncElement* expr,
                          const Array<MultiIndex>& mi,
                          Array<RCP<EvalVector> >& vec) const 
{
  static Array<string> coordNames;


  Tabs tabs;

  SUNDANCE_MSG1(verb(),
               tabs << "evaluating discrete func " << expr->toString() 
               << " with multiindices " << mi);

  if (coordNames.size() != 3)
    {
      coordNames.resize(3);
      coordNames[0] = "x";
      coordNames[1] = "y";
      coordNames[2] = "z";
    }



  std::string funcName = expr->name();
  
  TEST_FOR_EXCEPTION(!funcIdToFieldNumberMap_.containsKey(expr->fid().dofID()),
                     InternalError, "funcID " << expr->fid().dofID()
                     << " not found in TestEvalMediator funcID to field "
                     "map" << funcIdToFieldNumberMap_);

  int fieldIndex = funcIdToFieldNumberMap_.get(expr->fid().dofID());
  
  for (int i=0; i<mi.size(); i++)
    {
      if (mi[i].order()==0)
        {
          vec[i]->setString(funcName);
        }
      else
        {
          int dir = mi[i].firstOrderDirection();
          std::string deriv = "D[" + funcName + ", " + coordNames[dir] + "]";
          vec[i]->setString(deriv);
        }

      double * const xx = vec[i]->start();
      SUNDANCE_MSG3(verb(), "coeff=" << fields_[fieldIndex].coeff());
      xx[0] = fields_[fieldIndex].coeff() * evalDummyBasis(fieldIndex, mi[i]);
    }

  if (verb() > 0)
    {
      Out::os() << tabs << "results:" << std::endl;
      for (int i=0; i<mi.size(); i++)
        {
          Tabs tab1;
          Out::os() << tab1 << "mi=" << mi[i].toString() 
               << *vec[i] << std::endl;
        }
    }
}

double TestEvalMediator::evalDummyBasis(int m, const MultiIndex& mi) const
{
  TEST_FOR_EXCEPTION(mi.order() > 1, RuntimeError, 
                     "TestEvalMediator::evalDummyBasis found multiindex "
                     "order > 1. The bad multiindex was " << mi.toString());

  ADReal result = fields_[m].basis().evaluate(ADField::evalPoint());
  SUNDANCE_MSG3(verb(), "basis.value() " << result.value());
  SUNDANCE_MSG3(verb(), "basis.gradient() " << result.gradient());

  if (mi.order()==0)
    {
      return result.value();
    }
  else 
    {
      return result.gradient()[mi.firstOrderDirection()];
    }
}


