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

#include "SundanceTabs.hpp"
#include "SundanceRegionQuadCombo.hpp"

using namespace Sundance;
using namespace Teuchos;


RegionQuadCombo::RegionQuadCombo()
  : id_(-1), domain_(), quad_(), paramCurve_() , watch_()
{;}

RegionQuadCombo::RegionQuadCombo(const RCP<CellFilterStub>& domain,
  const RCP<QuadratureFamilyStub>& quad,
  const ParametrizedCurve& paramCurve ,
  const WatchFlag& watch)
  : id_(getID(domain, quad,watch)), domain_(domain), quad_(quad),
  paramCurve_(paramCurve) , watch_(watch)
{;}


int RegionQuadCombo::getID(const RCP<CellFilterStub>& domain,
  const RCP<QuadratureFamilyStub>& quad,
  const WatchFlag& watch)
{
  RegTriple p(domain, quad, watch);

  if (!domainAndQuadToIDMap().containsKey(p))
    {
      int id = topID();
      domainAndQuadToIDMap().put(p, id);
    }
  return domainAndQuadToIDMap().get(p);
}

string RegionQuadCombo::toString() const
{
  TeuchosOStringStream oss;
  Tabs tabs;

  oss << "Integration Region" << std::endl;
  {
    oss << tabs << "cell filter=" << domain_->description() << std::endl;
    oss << tabs << "quadrature rule=" << quad_->description() << std::endl;
    oss << tabs << "watchpoint=[" << watch().name() << "]" << std::endl;
  }
  return oss.str();
}

Sundance::Map<RegTriple, int>& RegionQuadCombo::domainAndQuadToIDMap()
{
  static Sundance::Map<RegTriple, int> rtn = Sundance::Map<RegTriple, int>();
  return rtn;
}





