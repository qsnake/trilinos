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

#include "SundanceMeshBase.hpp"
#include "SundanceMesh.hpp"
#include "SundanceExceptions.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceTabs.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using namespace Sundance;


MeshBase::MeshBase(int dim, const MPIComm& comm, 
  const MeshEntityOrder& order) 
  : dim_(dim), 
    comm_(comm),
    order_(order),
    reorderer_(Mesh::defaultReorderer().createInstance(this)),
    validWeights_(true),
    specialWeights_(),
    curvePoints_Are_Valid_(),
    nrCurvesForIntegral_(0),
    curvePoints_() ,
    curveDerivative_(),
    curveNormal_(),
    curveID_to_ArrayIndex_()
{ specialWeights_.resize(dim_);
  curvePoints_.resize(0);
  curveDerivative_.resize(0);
  curveNormal_.resize(0);
}



Point MeshBase::centroid(int cellDim, int cellLID) const
{
  if (cellDim==0) return nodePosition(cellLID);
  int dummy;
  Point x = nodePosition(facetLID(cellDim, cellLID, 0, 0, dummy));
  int nf = numFacets(cellDim, cellLID, 0);
  for (int f=1; f<nf; f++) 
    x += nodePosition(facetLID(cellDim, cellLID, 0, f, dummy));
  return x / ((double) nf);
}

void MeshBase::outwardNormals(
  const Array<int>& cellLIDs,
  Array<Point>& outwardNormals
  ) const 
{
  int D = spatialDim();
  outwardNormals.resize(cellLIDs.size());
  for (int c=0; c<cellLIDs.size(); c++)
  {
    int f=-1;
    TEST_FOR_EXCEPTION(numMaxCofacets(D-1, cellLIDs[c]) > 1, 
      RuntimeError,
      "cell #" << cellLIDs[c] << " is not a boundary cell");
    int maxLID = maxCofacetLID(D-1, cellLIDs[c], 0, f);
    Point cInterior = centroid(D, maxLID);
    Point cBdry = centroid(D-1, cellLIDs[c]);
    Point q = cBdry - cInterior;
    Point s;
    if (D==1) 
    {
      s = Point(1.0);
    }
    else if (D==2)
    {
      Point A = nodePosition(facetLID(D-1, cellLIDs[c], 0, 0, f));
      Point B = nodePosition(facetLID(D-1, cellLIDs[c], 0, 1, f));
      Point t = B - A;
      s = Point(-t[1], t[0]);
    }
    else 
    {
      Point A = nodePosition(facetLID(D-1, cellLIDs[c], 0, 0, f));
      Point B = nodePosition(facetLID(D-1, cellLIDs[c], 0, 1, f));
      Point C = nodePosition(facetLID(D-1, cellLIDs[c], 0, 2, f));
      s = cross(B-A, C-A);
    }
    if (q * s > 0.0)
    {
      outwardNormals[c] = s/::sqrt(s*s);
    }
    else
    {
      outwardNormals[c] = -s/::sqrt(s*s);
    }
  }
}


void  MeshBase::tangentsToEdges(
  const Array<int>& cellLIDs,
  Array<Point>& tangentVectors
  ) const 
{
  TEST_FOR_EXCEPT(spatialDim() <= 1);

  tangentVectors.resize(cellLIDs.size());

  for (int c=0; c<cellLIDs.size(); c++)
  {
    int fOrient=1;
    Point A = nodePosition(facetLID(1, cellLIDs[c], 0, 0, fOrient));
    Point B = nodePosition(facetLID(1, cellLIDs[c], 0, 1, fOrient));
    Point t = B - A;
    tangentVectors[c] = t/(sqrt(t*t));
  }
}




void MeshBase::getFacetArray(int cellDim, int cellLID, int facetDim, 
  Array<int>& facetLIDs,
  Array<int>& facetOrientations) const
{
  int nf = numFacets(cellDim, cellLID, facetDim);
  facetLIDs.resize(nf);
  facetOrientations.resize(nf);
  for (int f=0; f<nf; f++) 
  {
    facetLIDs[f] = facetLID(cellDim, cellLID, facetDim, f, 
      facetOrientations[f]);
  }
}


void MeshBase::getLabels(int cellDim, const Array<int>& cellLID, 
  Array<int>& labels) const
{
  labels.resize(cellLID.size());
  for (int i=0; i<cellLID.size(); i++) labels[i] = label(cellDim, cellLID[i]);
}


// ===================== storing curve intersection/quadrature points ======================

bool MeshBase::hasCurvePoints(int maxCellLID , int curveID) const {
   	if (curvePoints_[mapCurveID_to_Index(curveID)].containsKey(maxCellLID))
		return ( curvePoints_[mapCurveID_to_Index(curveID)].get(maxCellLID).size() > 0 );
   	else
   		return false;
}

void MeshBase::setCurvePoints(int maxCellLID, int curveID ,
		Array<Point>& points , Array<Point>& derivs , Array<Point>& normals) const {

	  Tabs tabs;
	  int verbo = 0;
	  SUNDANCE_MSG3(verbo, tabs << "MeshBase::setCurvePoints , nr:" << nrCurvesForIntegral_);
   	  curvePoints_[mapCurveID_to_Index(curveID)].put( maxCellLID , points );
   	  curveDerivative_[mapCurveID_to_Index(curveID)].put( maxCellLID , derivs );
   	  curveNormal_[mapCurveID_to_Index(curveID)].put( maxCellLID , normals );

}

void MeshBase::getCurvePoints(int maxCellLID, int curveID ,
		Array<Point>& points , Array<Point>& derivs , Array<Point>& normals) const {

	  Tabs tabs;
	  int verbo = 0;
	  SUNDANCE_MSG3(verbo, tabs << "MeshBase::getCurvePoints , nr:" << nrCurvesForIntegral_);
   	  points = curvePoints_[mapCurveID_to_Index(curveID)].get( maxCellLID );
   	  derivs = curveDerivative_[mapCurveID_to_Index(curveID)].get( maxCellLID );
   	  normals = curveNormal_[mapCurveID_to_Index(curveID)].get( maxCellLID );

}

int MeshBase::mapCurveID_to_Index(int curveID) const {

	 Tabs tabs;
	 int verbo = 0;

	 SUNDANCE_MSG3(verbo, tabs << "MeshBase::mapCurveID_to_Index curveID:" << curveID);
     if (curveID_to_ArrayIndex_.containsKey(curveID)){
    	 SUNDANCE_MSG3(verbo, tabs << "MeshBase::mapCurveID_to_Index value found for curveID:" << curveID << " ret:" << curveID_to_ArrayIndex_.get(curveID));
       	 return curveID_to_ArrayIndex_.get(curveID);
     } else {
    	 SUNDANCE_MSG3(verbo, tabs << "MeshBase::mapCurveID_to_Index create new :" << nrCurvesForIntegral_);
       	 curveID_to_ArrayIndex_.put( curveID , nrCurvesForIntegral_ );
       	 SUNDANCE_MSG3(verbo, tabs << "MeshBase::mapCurveID_to_Index , increment ");
       	 nrCurvesForIntegral_++;
    	 curvePoints_.resize(nrCurvesForIntegral_);
    	 curveDerivative_.resize(nrCurvesForIntegral_);
    	 curveNormal_.resize(nrCurvesForIntegral_);
    	 SUNDANCE_MSG3(verbo, tabs << "MeshBase::mapCurveID_to_Index create new :" << nrCurvesForIntegral_);
       	 return nrCurvesForIntegral_-1;
     }

}
