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

#include "SundanceQuadratureFamilyBase.hpp"


using namespace Sundance;
using namespace Teuchos;

void QuadratureFamilyBase::getPoints(const CellType& cellType, 
                                     Array<Point>& quadPoints,
                                     Array<double>& quadWeights) const
{
  switch(cellType)
    {
    case PointCell:
      quadPoints = tuple(Point());
      quadWeights = tuple(1.0);
      break;
    case LineCell:
      getLineRule(quadPoints, quadWeights);
      break;
    case TriangleCell:
      getTriangleRule(quadPoints, quadWeights);
      break;
    case QuadCell:
      getQuadRule(quadPoints, quadWeights);
      break;
    case TetCell:
      getTetRule(quadPoints, quadWeights);
      break;
    case BrickCell:
      getBrickRule(quadPoints, quadWeights);
      break;
    default:
#ifndef TRILINOS_7
      SUNDANCE_ERROR("cell type " << cellType << " not handled in "
                     "QuadratureFamilyBase::getPoints()");
#else
      SUNDANCE_ERROR7("cell type " << cellType << " not handled in "
                     "QuadratureFamilyBase::getPoints()");
#endif
    }
}

void QuadratureFamilyBase::getAdaptedWeights(const CellType& cellType ,
									 int cellDim,
	                                 int celLID ,
	                	             int facetIndex ,
	                                 const Mesh& mesh ,
	                                 const ParametrizedCurve& globalCurve ,
	                                 Array<Point>& quadPoints ,
	                                 Array<double>& quadWeights ,
	                                 bool &isCut) const {
#ifndef TRILINOS_7
  SUNDANCE_ERROR("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for " << toXML());
#else
  SUNDANCE_ERROR7("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for " << toXML());
#endif
}

void QuadratureFamilyBase::getLineRule(Array<Point>& /* quadPoints */,
                                       Array<double>& /* quadWeights */) const 
{
#ifndef TRILINOS_7
  SUNDANCE_ERROR("Line rule not available for " << toXML());
#else
  SUNDANCE_ERROR7("Line rule not available for " << toXML());
#endif
}

void QuadratureFamilyBase::getTriangleRule(Array<Point>& /* quadPoints */,
                                           Array<double>& /* quadWeights */) const 
{
#ifndef TRILINOS_7
  SUNDANCE_ERROR("Triangle rule not available for " << toXML());
#else
  SUNDANCE_ERROR7("Triangle rule not available for " << toXML());
#endif
}

void QuadratureFamilyBase::getQuadRule(Array<Point>& /* quadPoints */,
                                           Array<double>& /* quadWeights */) const 
{
#ifndef TRILINOS_7
  SUNDANCE_ERROR("Quad cell rule not available for " << toXML());
#else
  SUNDANCE_ERROR7("Quad cell rule not available for " << toXML());
#endif
}

void QuadratureFamilyBase::getTetRule(Array<Point>& /* quadPoints */,
                                           Array<double>& /* quadWeights */) const 
{
#ifndef TRILINOS_7
  SUNDANCE_ERROR("Tet cell rule not available for " << toXML());
#else
  SUNDANCE_ERROR7("Tet cell rule not available for " << toXML());
#endif
}

void QuadratureFamilyBase::getBrickRule(Array<Point>& /* quadPoints */,
                                           Array<double>& /* quadWeights */) const 
{
#ifndef TRILINOS_7
  SUNDANCE_ERROR("Brick rule not available for " << toXML());
#else
  SUNDANCE_ERROR7("Brick rule not available for " << toXML());
#endif
}
