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

#include "SundanceGaussianQuadratureType.hpp"
#include "SundanceGaussianQuadrature.hpp"
#include "SundanceGauss1D.hpp"
#include "SundanceTriangleQuadrature.hpp"
#include "SundanceQuadQuadrature.hpp"
#include "SundanceTetQuadrature.hpp"
#include "SundanceBrickQuadrature.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

GaussianQuadratureType::GaussianQuadratureType()
  : QuadratureTypeBase()
{
  
}

XMLObject GaussianQuadratureType::toXML() const 
{
  XMLObject rtn("GaussianQuadratureType");
  return rtn;
}


bool GaussianQuadratureType::supportsCellType(const CellType& cellType) const
{
  switch(cellType)
    {
    case PointCell:
    case LineCell:
    case TriangleCell:
    case QuadCell:
    case TetCell:
    case BrickCell:
      return true;
    default:
      return false;
    }
}

int GaussianQuadratureType::maxOrder(const CellType& cellType) const
{
  switch(cellType)
    {
    case TetCell:
      return TetQuadrature::maxOrder();
    default:
      return -1;
    }
}


bool GaussianQuadratureType::hasLimitedOrder(const CellType& cellType) const
{
  switch(cellType)
    {
    case TetCell:
      return true;
    default:
      return false;
    }
}

bool GaussianQuadratureType::supports(const CellType& cellType, int order) const
{
  if (order <= 0) return false;

  switch(cellType)
    {
    case PointCell:
      return true;
    case LineCell:
      return true;
    case TriangleCell:
      return true;
    case QuadCell:
      return true;
    case TetCell:
      return TetQuadrature::supportsOrder(order);
    case BrickCell:
       return true;
    default:
      return false;
    }
}

QuadratureFamily GaussianQuadratureType::createQuadFamily(int order) const
{
  return new GaussianQuadrature(order);
}
