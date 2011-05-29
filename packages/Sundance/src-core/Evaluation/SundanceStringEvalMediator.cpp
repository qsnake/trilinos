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


#include "SundanceStringEvalMediator.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceCellDiameterExpr.hpp"
#include "SundanceCellVectorExpr.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceExceptions.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace std;



StringEvalMediator::StringEvalMediator()
  : AbstractEvalMediator() 
{}

void StringEvalMediator::evalCoordExpr(const CoordExpr* expr,
                                       RCP<EvalVector>& vec) const
{
  SUNDANCE_MSG1(verb(), 
    "evaluating coord expr " << expr->toString());
  
  vec->resize(1);
  vec->start()[0]=0.0;
  vec->setString(expr->name());
}

void StringEvalMediator::evalCellDiameterExpr(const CellDiameterExpr* expr,
                                              RCP<EvalVector>& vec) const
{
  SUNDANCE_MSG1(verb(), 
    "evaluating cell diameter expr " << expr->toXML().toString());
  
  vec->resize(1);
  vec->start()[0]=0.0;
  vec->setString(expr->name());
}

void StringEvalMediator::evalCellVectorExpr(const CellVectorExpr* expr,
                                              RCP<EvalVector>& vec) const
{
  SUNDANCE_MSG1(verb(), "evaluating cell vector expr " << expr->toXML().toString());
  
  vec->resize(1);
  vec->start()[0]=0.0;
  vec->setString(expr->name());
}

void StringEvalMediator
::evalDiscreteFuncElement(const DiscreteFuncElement* expr,
                          const Array<MultiIndex>& mi,
                          Array<RCP<EvalVector> >& vec) const 
{
  static Array<string> coordNames;

  if (coordNames.size() != 3)
    {
      coordNames.resize(3);
      coordNames[0] = "x";
      coordNames[1] = "y";
      coordNames[2] = "z";
    }

  std::string funcName = expr->name();
  
  for (int i=0; i<mi.size(); i++)
    {
      vec[i]->resize(1);
      vec[i]->start()[0]=0.0;
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
    }
}
