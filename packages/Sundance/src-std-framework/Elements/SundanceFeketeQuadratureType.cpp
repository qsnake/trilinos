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

#include "SundanceFeketeQuadratureType.hpp"
#include "SundanceFeketeQuadrature.hpp"
#include "SundanceGaussLobatto1D.hpp"
#include "SundanceFeketeTriangleQuadrature.hpp"

using namespace Sundance;
using namespace Teuchos;

FeketeQuadratureType::FeketeQuadratureType() :
	QuadratureTypeBase()
{

}

XMLObject FeketeQuadratureType::toXML() const
{
	XMLObject rtn("FeketeQuadratureType");
	return rtn;
}

bool FeketeQuadratureType::supportsCellType(const CellType& cellType) const
{
	switch (cellType)
	{
	case PointCell:
	case LineCell:
	case TriangleCell:
	case QuadCell:
	case BrickCell:
		return true;
	default:
		return false;
	}
}

int FeketeQuadratureType::maxOrder(const CellType& cellType) const
{
	switch (cellType)
	{
	case TriangleCell:
		return FeketeTriangleQuadrature::maxOrder();
	default:
		return -1;
	}
}

bool FeketeQuadratureType::hasLimitedOrder(const CellType& cellType) const
{
	switch (cellType)
	{
	case TriangleCell:
		return true;
	default:
		return false;
	}
}

bool FeketeQuadratureType::supports(const CellType& cellType, int order) const
{
	if (order <= 0)
		return false;

	switch (cellType)
	{
	case PointCell:
		return true;
	case LineCell:
		return true;
	case QuadCell:
		return true;
	case BrickCell:
		return true;
	case TriangleCell:
		return FeketeTriangleQuadrature::supportsOrder(order);
	default:
		return false;
	}
}

QuadratureFamily FeketeQuadratureType::createQuadFamily(int order) const
{
	return new FeketeQuadrature(order);
}
