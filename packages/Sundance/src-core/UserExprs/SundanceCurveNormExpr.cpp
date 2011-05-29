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

#include "SundanceCurveNormExpr.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceOut.hpp"
#include "SundanceObjectWithVerbosity.hpp"

using namespace Sundance;
using namespace Teuchos;


CurveNormExpr::CurveNormExpr(int dir, const std::string& name)
  : EvaluatableExpr(),
    dir_(dir),
    name_(coordName(dir, name))
{}

bool CurveNormExpr::lessThan(const ScalarExpr* other) const
{
  const CurveNormExpr* c = dynamic_cast<const CurveNormExpr*>(other);
  TEST_FOR_EXCEPTION(c==0, InternalError, "cast should never fail at this point");
  return dir() < c->dir();
}

XMLObject CurveNormExpr::toXML() const
{
  XMLObject rtn("CurveNormExpr");
  rtn.addAttribute("dir", Teuchos::toString(dir_));
  rtn.addAttribute("name", name());
  return rtn;
}

string CurveNormExpr::coordName(int dir, const std::string& name)
{
  if (name.length() > 0) return name;
  switch(dir)
    {
    case 0:
      return "nx";
    case 1:
      return "ny";
    case 2:
      return "nz";
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError,
                         "CurveNormExpr::coordName direction out of range [0,2]");
      return "error";
    }
}


Set<MultipleDeriv> 
CurveNormExpr::internalFindW(int order, const EvalContext& context) const
{
  Tabs tab0;
  SUNDANCE_VERB_HIGH(tab0 << "CurveNormExpr::internalFindW() for " << toString());
  Set<MultipleDeriv> rtn;

  if (order==0) rtn.put(MultipleDeriv());

  if (order==1) 
    {
	  /* todo: in case of first derivative currently we return nothing,
       * considering the normal as constant, which mathematically is not
       * correct.
       * However calculating the derivative of the normal vector in one
       * direction is not trivial, remains to be done later */
    }

  SUNDANCE_VERB_HIGH(tab0 << "W[" << order << "]=" << rtn);
  return rtn;
}

/*
Set<MultipleDeriv> 
CurveNormEvaluator::internalFindV(int order, const EvalContext& context) const
{
  Tabs tab0;
  SUNDANCE_VERB_HIGH(tab0 << "CurveNormEvaluator::internalFindV() for " << toString());
  Set<MultipleDeriv> rtn;

  if (order==0) rtn.put(MultipleDeriv());

  SUNDANCE_VERB_HIGH(tab0 << "V[" << order << "]=" << rtn);
  return rtn;
}
*/

/*
Set<MultipleDeriv> 
CurveNormEvaluator::internalFindC(int order, const EvalContext& context) const
{
  Tabs tab0;
  SUNDANCE_VERB_HIGH(tab0 << "CurveNormEvaluator::internalFindC() for " << toString());
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
}*/






