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

#include "SundanceSumOfBCs.hpp"
#include "SundanceTabs.hpp"

using namespace Sundance;
using namespace Teuchos;
using std::endl;

SumOfBCs::SumOfBCs(const RCP<CellFilterStub>& region,
  const Expr& expr,
  const RCP<QuadratureFamilyStub>& quad,
  const ParametrizedCurve& curve,
  const WatchFlag& watch)
  : SumOfIntegrals(region, expr, quad, curve, watch)
{;}

std::ostream& SumOfBCs::toText(std::ostream& os, bool paren) const
{
  os << "Sum of BCs[" << std::endl;
  for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap().begin(); i!=rqcToExprMap().end(); i++)
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


XMLObject SumOfBCs::toXML() const 
{
  XMLObject rtn("SumOfBCs");
  for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap().begin(); i!=rqcToExprMap().end(); i++)
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
