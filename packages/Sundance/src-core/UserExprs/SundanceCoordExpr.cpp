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

#include "SundanceCoordExpr.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceOut.hpp"
#include "SundanceObjectWithVerbosity.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;


CoordExpr::CoordExpr(int dir, const std::string& name)
  : EvaluatableExpr(),
    dir_(dir),
    name_(coordName(dir, name))
{}

bool CoordExpr::lessThan(const ScalarExpr* other) const
{
  const CoordExpr* c = dynamic_cast<const CoordExpr*>(other);
  TEST_FOR_EXCEPTION(c==0, InternalError, "cast should never fail at this point");
  return dir() < c->dir();
}

XMLObject CoordExpr::toXML() const 
{
  XMLObject rtn("CoordExpr");
  rtn.addAttribute("dir", Teuchos::toString(dir_));
  rtn.addAttribute("name", name());
  return rtn;
}

string CoordExpr::coordName(int dir, const std::string& name)
{
  if (name.length() > 0) return name;
  switch(dir)
    {
    case 0:
      return "x";
    case 1:
      return "y";
    case 2:
      return "z";
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError,
                         "CoordExpr::coordName direction out of range [0,2]");
      return "error";
    }
}


Set<MultipleDeriv> 
CoordExpr::internalFindW(int order, const EvalContext& context) const
{
  Tabs tab0;
  SUNDANCE_VERB_HIGH(tab0 << "CoordExpr::internalFindW() for " << toString());  
  Set<MultipleDeriv> rtn;

  if (order==0) rtn.put(MultipleDeriv());

  if (order==1) 
    {
      Deriv x = coordDeriv(dir_);
      MultipleDeriv md;
      md.put(x);
      rtn.put(md);
    }

  SUNDANCE_VERB_HIGH(tab0 << "W[" << order << "]=" << rtn);
  return rtn;
}


Set<MultipleDeriv> 
CoordExpr::internalFindV(int order, const EvalContext& context) const
{
  Tabs tab0;
  SUNDANCE_VERB_HIGH(tab0 << "CoordExpr::internalFindV() for " << toString());  
  Set<MultipleDeriv> rtn;

  if (order==0) rtn.put(MultipleDeriv());

  SUNDANCE_VERB_HIGH(tab0 << "V[" << order << "]=" << rtn);
  return rtn;
}


Set<MultipleDeriv> 
CoordExpr::internalFindC(int order, const EvalContext& context) const
{
  Tabs tab0;
  SUNDANCE_VERB_HIGH(tab0 << "CoordExpr::internalFindC() for " << toString());  
  Set<MultipleDeriv> rtn;

  if (order==1) 
    {
      Deriv x = coordDeriv(dir_);
      MultipleDeriv md;
      md.put(x);
      rtn.put(md);
    }
  SUNDANCE_VERB_HIGH(tab0 << "C[" << order << "]=" << rtn);
  return rtn;
}






