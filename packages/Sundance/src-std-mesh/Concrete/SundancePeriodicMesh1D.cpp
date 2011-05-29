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

#include "SundancePeriodicMesh1D.hpp"

#include "SundanceCellJacobianBatch.hpp"
#include "SundanceMaximalCofacetBatch.hpp"
#include "SundanceDebug.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_MPIContainerComm.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceObjectWithVerbosity.hpp"

using namespace Sundance;
using namespace Teuchos;
using namespace std;


PeriodicMesh1D::PeriodicMesh1D(double xMin, double xMax, int numElems)
  : MeshBase(1, MPIComm::self(), ExodusMeshOrder),
    numElems_(numElems),
    xMin_(xMin),
    xMax_(xMax),
    x_(numElems+1),
    verts_(numElems),
    labels_(2)
{
  labels_[0].resize(numElems_);
  labels_[1].resize(numElems_);

  for (int i=0; i<numElems_; i++)
  {
    labels_[0][i] = 0;
    labels_[1][i] = 1;
    verts_[i].resize(2);
    verts_[i][0] = i;
    verts_[i][1] = (i+1) % numElems;
  }

  for (int i=0; i<=numElems_; i++)
  {
    x_[i] = Point(xMin_ + ((double) i)/((double) numElems_)*(xMax_-xMin_));
  }
}

int PeriodicMesh1D::numCells(int cellDim) const
{
  switch(cellDim)
  {
    case 0 :
      return numElems_;
    case 1:
      return numElems_;
    default:
      TEST_FOR_EXCEPT(true);
  }
  return -1; // -Wall
}


Point PeriodicMesh1D::nodePosition(int i) const 
{
  return x_[i];
}


const double* PeriodicMesh1D::nodePositionView(int i) const 
{
  return &(x_[i][0]);
}

void PeriodicMesh1D::getJacobians(int cellDim, const Array<int>& cellLID,
    CellJacobianBatch& jBatch) const
{
  TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), InternalError,
    "cellDim=" << cellDim 
    << " is not in expected range [0, " << spatialDim()
    << "]");

  jBatch.resize(cellLID.size(), spatialDim(), cellDim);

  int nCells = cellLID.size();

  if (cellDim==0)
  {
    for (int i=0; i<nCells; i++)
    {
      double* detJ = jBatch.detJ(i);
      *detJ = 1.0;
    }
  }
  else
  {
    for (int i=0; i<nCells; i++)
    {
      int lid = cellLID[i];
      double* J = jBatch.jVals(i);
      double x0 = x_[lid][0];
      double x1 = x_[lid+1][0];
      J[0] = fabs(x1-x0);
    }
  }
}


void PeriodicMesh1D::getCellDiameters(int cellDim, const Array<int>& cellLID,
  Array<double>& cellDiameters) const
{
  cellDiameters.resize(cellLID.size());
  
  TEST_FOR_EXCEPT(cellDim < 0 || cellDim > 1);

  if (cellDim==1)
  {
    for (int i=0; i<cellLID.size(); i++)
    {
      int c = cellLID[i];
      double h = x_[c+1][0]-x_[c][0];
      cellDiameters[i] = h;
    }
  }
  else
  {
    for (int i=0; i<cellLID.size(); i++)
    {
      cellDiameters[i] = 1.0;
    }
  }
}

void PeriodicMesh1D::pushForward(int cellDim, const Array<int>& cellLID,
    const Array<Point>& refQuadPts,
    Array<Point>& physQuadPts) const
{
  TEST_FOR_EXCEPT(cellDim < 0 || cellDim > 1);

  if (cellDim==1)
  {
    if (physQuadPts.size() > 0) physQuadPts.resize(0);
    physQuadPts.reserve(refQuadPts.size() * cellLID.size());
    
    for (int i=0; i<cellLID.size(); i++)
    {
      int c = cellLID[i];
      Point h = x_[c+1]-x_[c];
      for (int q=0; q<refQuadPts.size(); q++)
      {
        physQuadPts.append(Point(x_[c] + refQuadPts[q][0] * h));
      }
    }
  }
  else
  {
    for (int i=0; i<cellLID.size(); i++)
    {
      physQuadPts.append(x_[cellLID[i]]);
    }
  }
}

int PeriodicMesh1D::numFacets(int cellDim, int cellLID,
    int facetDim) const
{
  TEST_FOR_EXCEPT(cellLID < 0 || cellLID >= numElems_);
  if (cellDim == 1 && facetDim==0) return 2;
  return 0;
}


    
int PeriodicMesh1D::facetLID(int cellDim, int cellLID,
  int facetDim, int facetIndex,
  int& facetOrientation) const
{
  TEST_FOR_EXCEPT(cellLID < 0 || cellLID >= numElems_);

  TEST_FOR_EXCEPT(cellDim != 1);
  TEST_FOR_EXCEPT(facetDim != 0);
  TEST_FOR_EXCEPT(facetIndex < 0);
  TEST_FOR_EXCEPT(facetIndex > 1);

  return verts_[cellLID][facetIndex];
}

void PeriodicMesh1D::getFacetLIDs(int cellDim,
    const Array<int>& cellLID,
    int facetDim,
    Array<int>& facetLID,
    Array<int>& facetSign) const
{
  facetLID.resize(2*cellLID.size());
  facetSign.resize(2*cellLID.size());

  for (int i=0; i<cellLID.size(); i++) 
  {
    facetLID[2*i] = this->facetLID(cellDim, cellLID[i], facetDim, 0, facetSign[2*i]);
    facetLID[2*i+1] = this->facetLID(cellDim, cellLID[i], facetDim, 1, facetSign[2*i]);
  }
}


const int* PeriodicMesh1D::elemZeroFacetView(int cellLID) const
{
  return &(verts_[cellLID][0]);
}


int PeriodicMesh1D::numMaxCofacets(int cellDim, int cellLID) const
{
  TEST_FOR_EXCEPT(cellDim != 0);
  return 2;
}

int PeriodicMesh1D::maxCofacetLID(int cellDim, int cellLID,
    int cofacetIndex,
    int& facetIndex) const
{
  TEST_FOR_EXCEPT(cellDim != 0);
  if (cofacetIndex==0) return (cellLID-1) % numElems_;
  return cellLID;
}

void PeriodicMesh1D::getMaxCofacetLIDs(const Array<int>& cellLIDs,
    MaximalCofacetBatch& cofacets) const
{
  TEST_FOR_EXCEPT(true);
}

void PeriodicMesh1D::getCofacets(int cellDim, int cellLID,
    int cofacetDim, Array<int>& cofacetLIDs) const
{
  TEST_FOR_EXCEPT(cellDim != 0);
  TEST_FOR_EXCEPT(cofacetDim != 1);

  cofacetLIDs.resize(2);
  cofacetLIDs[0] = (cellLID-1) % numElems_;
  cofacetLIDs[1] = cellLID;
}

int PeriodicMesh1D::mapGIDToLID(int cellDim, int globalIndex) const
{
  return globalIndex;
}

bool PeriodicMesh1D::hasGID(int cellDim, int globalIndex) const
{
  return true;
}

int PeriodicMesh1D::mapLIDToGID(int cellDim, int localIndex) const
{
  return localIndex;
}

CellType PeriodicMesh1D::cellType(int cellDim) const
{
  if (cellDim==0) return PointCell;
  else if (cellDim==1) return LineCell;
  else return NullCell;
}

int PeriodicMesh1D::label(int cellDim, int cellLID) const
{
  return labels_[cellDim][cellLID];
}

void PeriodicMesh1D::getLabels(int cellDim, const Array<int>& cellLID,
  Array<int>& labels) const
{
  labels.resize(cellLID.size());
  const Array<int>& ld = labels_[cellDim];
  for (int i=0; i<cellLID.size(); i++)
  {
    labels[i] = ld[cellLID[i]];
  }
}

Set<int> PeriodicMesh1D::getAllLabelsForDimension(int cellDim) const
{
  Set<int> rtn;

  const Array<int>& ld = labels_[cellDim];
  for (int i=0; i<ld.size(); i++)
  {
    rtn.put(ld[i]);
  }
  
  return rtn;
}

void PeriodicMesh1D::setLabel(int cellDim, int cellLID, int label)
{
  labels_[cellDim][cellLID] = label;
}


void PeriodicMesh1D::getLIDsForLabel(int cellDim, int label, Array<int>& cellLIDs) const
{
  cellLIDs.resize(0);
  const Array<int>& ld = labels_[cellDim];
  for (int i=0; i<ld.size(); i++)
  {
    if (ld[i] == label) cellLIDs.append(i);
  }
}
