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

#include "SundanceClosedNewtonCotes.hpp"
#include "SundanceTriangleQuadrature.hpp"
#include "SundanceTetQuadrature.hpp"

using namespace Sundance;
using namespace Teuchos;

ClosedNewtonCotes::ClosedNewtonCotes(int order)
  : QuadratureFamilyBase(order)
{
  TEST_FOR_EXCEPTION(!(order==2 || order==3), RuntimeError, "order " 
    << order << " not supported by ClosedNewtonCotes");
}

XMLObject ClosedNewtonCotes::toXML() const 
{
  XMLObject rtn("ClosedNewtonCotes");
  rtn.addAttribute("order", Teuchos::toString(order()));
  return rtn;
}


void ClosedNewtonCotes::getLineRule(Array<Point>& quadPoints,
                                     Array<double>& quadWeights) const 
{
  Array<double> x;
  if (order()==2)
  {
    x = tuple(0.0, 1.0);
    quadWeights = tuple(0.5, 0.5);
  }
  else 
  {
    x = tuple(0.0, 0.5, 1.0);
    quadWeights = tuple(1.0/6.0, 2.0/3.0, 1.0/6.0);
  }
  quadPoints.resize(x.size());
    
  for (int i=0; i<x.length(); i++)
  {
    quadPoints[i] = Point(x[i]);
  }
}

void ClosedNewtonCotes::getTriangleRule(Array<Point>& quadPoints,
                                          Array<double>& quadWeights) const 
{
  if (order()==2)
  {
    quadPoints=tuple(Point(0.0, 0.0), Point(1.0, 0.0), Point(0.0, 1.0));
    quadWeights=tuple(1.0/3.0, 1.0/3.0, 1.0/3.0);
  }
  else
  {
    quadPoints = tuple(Point(0.0, 0.0), Point(0.5, 0.0), Point(1.0, 0.0),
      Point(0.5, 0.5), Point(0.0, 1.0), Point(0.0, 0.5),
      Point(1.0/3.0, 1.0/3.0));
    quadWeights = tuple(3.0/60.0, 8.0/60.0, 3.0/60.0,
      8.0/60.0, 3.0/60.0, 8.0/60.0,
      27.0/60.0);
  }
	for (int i=0; i<quadWeights.size(); i++)
		{
			quadWeights[i] *=  0.5;
		}
}

