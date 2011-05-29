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

#include "SundanceGaussianQuadrature.hpp"
#include "SundanceGauss1D.hpp"
#include "SundanceTriangleQuadrature.hpp"
#include "SundanceQuadQuadrature.hpp"
#include "SundanceTetQuadrature.hpp"
#include "SundanceBrickQuadrature.hpp"

using namespace Sundance;
using namespace Teuchos;

GaussianQuadrature::GaussianQuadrature(int order)
  : QuadratureFamilyBase(order)
{
  
}

XMLObject GaussianQuadrature::toXML() const 
{
  XMLObject rtn("GaussianQuadrature");
  rtn.addAttribute("order", Teuchos::toString(order()));
  return rtn;
}



void GaussianQuadrature::getLineRule(Array<Point>& quadPoints,
                                     Array<double>& quadWeights) const 
{
  int p = order() + 1;
  p = p + (p%2);
  int n = p/2;
			
  quadPoints.resize(n);
  quadWeights.resize(n);
			
  Gauss1D q1(n, 0.0, 1.0);
			
  for (int i=0; i<n; i++)
    {
      quadWeights[i] = q1.weights()[i];
      quadPoints[i] = Point(q1.nodes()[i]);
    }
}

void GaussianQuadrature::getTriangleRule(Array<Point>& quadPoints,
                                          Array<double>& quadWeights) const 
{
  Array<double> x;
  Array<double> y;
  Array<double> w;
			
  TriangleQuadrature::getPoints(order(), w, x, y);
  quadPoints.resize(w.length());
  quadWeights.resize(w.length());
  for (int i=0; i<w.length(); i++)
    {
      quadWeights[i] = 0.5*w[i];
      quadPoints[i] = Point(x[i], y[i]);
    }  
}

void GaussianQuadrature::getQuadRule(Array<Point>& quadPoints,
                                          Array<double>& quadWeights) const
{
  Array<double> x;
  Array<double> y;
  Array<double> w;

  QuadQuadrature::getPoints(order(), w, x, y);
  quadPoints.resize(w.length());
  quadWeights.resize(w.length());
  for (int i=0; i<w.length(); i++)
    {
      quadWeights[i] = w[i];
      quadPoints[i] = Point(x[i], y[i]);
    }
}

void GaussianQuadrature::getTetRule(Array<Point>& quadPoints,
                                    Array<double>& quadWeights) const 
{
  Array<double> x;
  Array<double> y;
  Array<double> z;
  Array<double> w;
			
  TetQuadrature::getPoints(order(), w, x, y, z);
  quadPoints.resize(w.length());
  quadWeights.resize(w.length());
  for (int i=0; i<w.length(); i++)
    {
      quadWeights[i] = w[i]/6.0;
      quadPoints[i] = Point(x[i], y[i], z[i]);
    }  
}

void GaussianQuadrature::getBrickRule(Array<Point>& quadPoints,
                                    Array<double>& quadWeights) const
{
  Array<double> x;
  Array<double> y;
  Array<double> z;
  Array<double> w;

  BrickQuadrature::getPoints(order(), w, x, y, z);
  quadPoints.resize(w.length());
  quadWeights.resize(w.length());
  for (int i=0; i<w.length(); i++)
    {
      quadWeights[i] = w[i];
      quadPoints[i] = Point(x[i], y[i], z[i]);
    }
}

void GaussianQuadrature::getAdaptedWeights(const CellType& cellType ,
									       int cellDim,
	                                       int celLID ,
	                	                   int facetIndex ,
	                                       const Mesh& mesh ,
	                                       const ParametrizedCurve& globalCurve ,
	                                       Array<Point>& quadPoints ,
	                                       Array<double>& quadWeights ,
	                                       bool &isCut) const {
	// First implementation ...

	// we say that the cell is always cut by the curve so these weights will be used
	isCut = true;
	Array<Point> tmpPoints = quadPoints;

    Array<int> celLIDs(1);
    celLIDs[0] = celLID;
    // transform the ref points to real coordinates
	mesh.pushForward( cellDim, celLIDs, quadPoints, tmpPoints );

	quadWeights.resize(tmpPoints.size());
    // simple weight calculation
	for (int i=0 ; i < tmpPoints.size() ; i++){
		quadWeights[i] = quadWeights[i] * globalCurve.integrationParameter(tmpPoints[i]);
	}

/*
	Array<Point> vertices;
    Array<int>  nodeLIDs;
    Array<int>  orient;
    int testcount = 1;

    // Get coordinates
    mesh.getFacetArray(cellDim, celLID, 0, nodeLIDs, orient);

    vertices.resize(nodeLIDs.size());
    for (int i=0; i<nodeLIDs.size(); i++) {
    	vertices[i] = mesh.nodePosition(nodeLIDs[i]);
    	SUNDANCE_MSG5(6, "Points [" << testcount << "/" << i << "] " << vertices[i]);
    	testcount++;
    }*/
}
