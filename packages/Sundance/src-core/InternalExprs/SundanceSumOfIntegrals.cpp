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

#include "SundanceSumOfIntegrals.hpp"
#include "SundanceSpatiallyConstantExpr.hpp"
#include "SundanceSpectralPreprocessor.hpp"
#include "SundanceTabs.hpp"

using namespace Sundance;
using namespace Teuchos;

SumOfIntegrals::SumOfIntegrals(const RCP<CellFilterStub>& region,
  const Expr& expr,
  const RCP<QuadratureFamilyStub>& quad,
  const WatchFlag& watch)
  : ScalarExpr(), rqcToExprMap_()
{
   addTerm(region, expr, quad, ParametrizedCurve::returnDummyCurve() , watch, 1);
}

SumOfIntegrals::SumOfIntegrals(const RCP<CellFilterStub>& region,
  const Expr& expr,
  const RCP<QuadratureFamilyStub>& quad,
  const ParametrizedCurve& curve,
  const WatchFlag& watch)  : ScalarExpr(), rqcToExprMap_()
{
  addTerm(region, expr, quad, curve , watch, 1);
}


Expr SumOfIntegrals::filterSpectral(const Expr& expr) const 
{
  return SpectralPreprocessor::projectSpectral(expr);
}



void SumOfIntegrals::addTerm(const RCP<CellFilterStub>& regionPtr,
  const Expr& expr,
  const RCP<QuadratureFamilyStub>& quadPtr, 
  const ParametrizedCurve& paramCurve,
  const WatchFlag& watch, int sign)
{
  Expr ex = filterSpectral(expr);

  RegionQuadCombo rqc(regionPtr, quadPtr, paramCurve ,watch);

  if (rqcToExprMap_.containsKey(rqc))
  {
    Expr e = rqcToExprMap_.get(rqc);
    rqcToExprMap_.put(rqc, e + sign*ex);
  }
  else
  {
    rqcToExprMap_.put(rqc, sign*ex);
  }
}


void SumOfIntegrals::merge(const SumOfIntegrals* other, int sign) 
{
  for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
         i=other->rqcToExprMap_.begin(); i!=other->rqcToExprMap_.end(); i++)
  {
    const RegionQuadCombo& rqc = i->first;
    const Expr& e = i->second;
    addTerm(rqc.domain(), e, rqc.quad(), rqc.paramCurve() , rqc.watch(), sign);
  }
}

void SumOfIntegrals::multiplyByConstant(const SpatiallyConstantExpr* expr) 
{
  double a = expr->value();
  Sundance::Map<RegionQuadCombo, Expr> newMap;
  for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    Expr e = i->second;
    newMap.put(i->first, a*e);
  }
  rqcToExprMap_ = newMap;
}

void SumOfIntegrals::changeSign()
{
  Sundance::Map<RegionQuadCombo, Expr> newMap;
  for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    Expr e = i->second;
    newMap.put(i->first, -e);
  }
  rqcToExprMap_ = newMap;
}

Set<int> SumOfIntegrals::funcsOnRegion(const OrderedHandle<CellFilterStub>& d, const Set<int>& funcSet) const 
{
  Set<int> rtn;
  for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    const RegionQuadCombo& rqc = i->first;
    if (OrderedHandle<CellFilterStub>(rqc.domain()) != d) continue;
    Expr e = i->second;
    e.ptr()->accumulateFuncSet(rtn, funcSet);
  }
  return rtn;
}


bool SumOfIntegrals::integralHasTestFunctions(const OrderedHandle<CellFilterStub>& d) const 
{
  for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    const RegionQuadCombo& rqc = i->first;
    if (OrderedHandle<CellFilterStub>(rqc.domain()) != d) continue;
    Expr e = i->second;
    if (e.hasTestFunctions()) return true;
  }
  return false;
}



RCP<CellFilterStub> SumOfIntegrals::nullRegion() const
{
  for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    const RegionQuadCombo& rqc = i->first;
    if (!rqc.domain()->isNullRegion())
    {
      return rqc.domain()->makeNullRegion();
    }
  }
  
  TEST_FOR_EXCEPTION(true, RuntimeError,
                     "SumOfIntegrals::nullRegion() called on a sum "
                     "of integrals with no non-null regions");

  return RCP<CellFilterStub>();
}

bool SumOfIntegrals::isIndependentOf(const Expr& u) const
{
  for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    Expr e = i->second;
    if (!e.isIndependentOf(u)) return false;
  }
  return true;
}

bool SumOfIntegrals::isLinearForm(const Expr& u) const
{
  for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    Expr e = i->second;
    if (!e.isLinearForm(u)) return false;
  }
  return true;
}

bool SumOfIntegrals::isQuadraticForm(const Expr& u) const
{
  for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    Expr e = i->second;
    if (!e.isQuadraticForm(u)) return false;
  }
  return true;
}


bool SumOfIntegrals::everyTermHasTestFunctions() const
{
  for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    Expr e = i->second;
    if (!e.everyTermHasTestFunctions()) return false;
  }
  return true;
}


bool SumOfIntegrals::isLinearInTests() const
{
  for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    Expr e = i->second;
    if (!e.isLinearInTests()) return false;
  }
  return true;
}

bool SumOfIntegrals::hasTestFunctions() const
{
  for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    Expr e = i->second;
    if (e.hasTestFunctions()) return true;
  }
  return false;
}



std::ostream& SumOfIntegrals::toText(std::ostream& os, bool paren) const
{
  os << "Sum of Integrals[" << std::endl;
  for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    const RegionQuadCombo& rqc = i->first;
    Expr e = i->second;
    os << "Integral[" << std::endl;
    os << "rqc=" << rqc.toString() << std::endl;
    os << "expr=" << e.toString() << std::endl;
    os << "]" << std::endl;
  }
  os << "]" << std::endl;

  return os;
}


XMLObject SumOfIntegrals::toXML() const 
{
  XMLObject rtn("SumOfIntegrals");
  for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    const RegionQuadCombo& rqc = i->first;
    Expr e = i->second;
    XMLObject child("Integral");
    rtn.addChild(child);
    child.addChild(rqc.quad()->toXML());
    child.addChild(rqc.domain()->toXML());
    child.addChild(rqc.watch().toXML());
    child.addChild(e.toXML());
  }

  return rtn;
}


bool SumOfIntegrals::lessThan(const ScalarExpr* other) const
{
  const SumOfIntegrals* f = dynamic_cast<const SumOfIntegrals*>(other);
  TEST_FOR_EXCEPTION(f==0, InternalError, "cast should never fail at this point");
  
  return rqcToExprMap_ < f->rqcToExprMap_;
}

bool SumOfIntegrals::hasWatchedTerm() const 
{
  for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    const RegionQuadCombo& rqc = i->first;
    if (rqc.watch().isActive()) return true;
  }
  return false;
}



int SumOfIntegrals::eqnSetSetupVerb() const 
{
  int rtn = 0;
  for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    const RegionQuadCombo& rqc = i->first;
    int v = rqc.watch().param("equation set setup");
    if (v > rtn) rtn = v;
  }
  return rtn;
}


