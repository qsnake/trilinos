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
#include "SundanceClosedNewtonCotesType.hpp"
#include "SundanceTriangleQuadrature.hpp"
#include "SundanceTetQuadrature.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

ClosedNewtonCotesType::ClosedNewtonCotesType()
  : QuadratureTypeBase()
{;}

XMLObject ClosedNewtonCotesType::toXML() const 
{
  XMLObject rtn("ClosedNewtonCotes");
  return rtn;
}

bool ClosedNewtonCotesType::supportsCellType(const CellType& cellType) const
{
  switch(cellType)
    {
    case PointCell:
    case LineCell:
    case TriangleCell:
      return true;
    default:
      return false;
    }
}

int ClosedNewtonCotesType::maxOrder(const CellType& cellType) const
{
  return 3;
}

bool ClosedNewtonCotesType::supports(const CellType& cellType, int order) const
{
  if (order <2 || order > 3) return false;
  switch(cellType)
    {
    case PointCell:
    case LineCell:
    case TriangleCell:
      return true;
    default:
      return false;
    }
}



QuadratureFamily ClosedNewtonCotesType::createQuadFamily(int order) const
{
  return new ClosedNewtonCotes(order);
}
