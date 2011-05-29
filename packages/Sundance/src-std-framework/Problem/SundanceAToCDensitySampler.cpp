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

#include "SundanceAToCDensitySampler.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceGeomUtils.hpp"
#include "SundanceLagrange.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include <queue>

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using namespace TSFExtended;


static Time& densitySamplingTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("density sampling"); 
  return *rtn;
}

AToCDensitySampler::AToCDensitySampler(const AToCPointLocator& locator,
                                       const VectorType<double>& vecType)
  : discSpace_(locator.mesh(), new Lagrange(0), locator.subdomain(), vecType),
    dim_(locator.mesh().spatialDim()),
    mesh_(locator.mesh()),
    elemToVecIndexMap_(),
    elemWeights_(new DiscreteFunction(discSpace_, 0.0)),
    elemWeightVec_(),
    locator_(locator),
    isAxisymmetric_(false),
    origin_(),
    axis_()
{
  init();
}


AToCDensitySampler::AToCDensitySampler(const AToCPointLocator& locator,
                                       const std::vector<double>& origin,
                                       const std::vector<double>& rotationalAxis,
                                       const VectorType<double>& vecType)
  : discSpace_(locator.mesh(), new Lagrange(0), locator.subdomain(), vecType),
    dim_(locator.mesh().spatialDim()),
    mesh_(locator.mesh()),
    elemToVecIndexMap_(),
    elemWeights_(new DiscreteFunction(discSpace_, 0.0)),
    elemWeightVec_(),
    locator_(locator),
    isAxisymmetric_(true),
    origin_(vec2point(origin)),
    axis_(normPoint(vec2point(rotationalAxis)))
{
  init();
}



void AToCDensitySampler::init()
{
  const CellFilter& domain = discSpace_.cellFilters(0);

  elemWeightVec_ = DiscreteFunction::discFunc(elemWeights_)->getVector();

  elemToVecIndexMap_ = rcp(new Array<int>(mesh_.numCells(dim_), -1));

  Array<int>& a = *elemToVecIndexMap_;

  CellSet cells = domain.getCells(mesh_);

  Array<int> cellLID;
  cellLID.reserve(mesh_.numCells(dim_));

  for (CellIterator i=cells.begin(); i!=cells.end(); i++)
    {
      cellLID.append(*i);
    }

  const RCP<DOFMapBase>& dofMap = discSpace_.map();

  Set<int> funcs = makeSet(0);
  Array<Array<int> > dofs;
  Array<int> nNodes;
  dofMap->getDOFsForCellBatch(dim_, cellLID, funcs, dofs, nNodes,0);
  
  const Array<int>& dofs0 = dofs[0];
  for (int c=0; c<cellLID.size(); c++)
    {
      int vecIndex = dofs0[c];
      int lid = cellLID[c];
      a[lid] = vecIndex;
      double vol = volume(mesh_, dim_, lid);
      if (isAxisymmetric_)
        {
          Point xCell = mesh_.centroid(dim_, lid) - origin_;
          double dPerp = ::sqrt(xCell*xCell - (xCell*axis_)*(xCell*axis_));
          vol = vol * dPerp;
        }
      elemWeightVec_.setElement(vecIndex, vol);
    }
}

Point AToCDensitySampler::vec2point(const std::vector<double>& x) const
{
  if (x.size()==1) return Point(x[0]);
  else if (x.size()==2U) return Point(x[0], x[1]);
  else if (x.size()==3U) return Point(x[0], x[1], x[2]);
  TEST_FOR_EXCEPT(x.size() < 1 || x.size() > 3U);
  return Point();
}

Point AToCDensitySampler::normPoint(const Point& x) const
{
  return x/sqrt(x*x);
}


Expr AToCDensitySampler::sample(const std::vector<double>& positions,
                                const double& particleWeight) const 
{
  TimeMonitor timer(densitySamplingTimer());

  TEST_FOR_EXCEPTION(positions.size() % dim_ != 0, RuntimeError,
                     "vector of coordinates should by an integer multiple "
                     "of the spatial dimension");

  Expr rtn = new DiscreteFunction(discSpace_, 0.0);
  Vector<double> density = DiscreteFunction::discFunc(rtn)->getVector();

  int nPts = positions.size() / dim_;

  for (int i=0; i<nPts; i++)
    {
      const double* x = &(positions[dim_*i]);

      int guess = locator_.guessCell(x);

      TEST_FOR_EXCEPTION(guess < 0, RuntimeError, "particle #" << i << " position="
                         << AToCPointLocator::makePoint(dim_, x) 
                         << " is out of search grid");

      int cellLID = locator_.findEnclosingCell(guess, x);

      int vecIndex = (*elemToVecIndexMap_)[cellLID];
      double vol = elemWeightVec_.getElement(vecIndex);
      density.addToElement(vecIndex, particleWeight/vol);
    }

  return rtn;
}


Expr AToCDensitySampler::resetCounts() const 
{
  Expr rtn = new DiscreteFunction(discSpace_, 0.0);

  return rtn;
}

void AToCDensitySampler::addToCounts(const std::vector<double>& positions,
                                     const double& particleWeight,
                                     Expr density) const 
{
  TimeMonitor timer(densitySamplingTimer());

  TEST_FOR_EXCEPTION(positions.size() % dim_ != 0, RuntimeError,
                     "vector of coordinates should by an integer multiple "
                     "of the spatial dimension");

  Vector<double> vec = DiscreteFunction::discFunc(density)->getVector();

  int nPts = positions.size() / dim_;

  for (int i=0; i<nPts; i++)
    {
      const double* x = &(positions[dim_*i]);

      int guess = locator_.guessCell(x);

      TEST_FOR_EXCEPTION(guess < 0, RuntimeError, "particle #" << i << " position="
                         << AToCPointLocator::makePoint(dim_, x) 
                         << " is out of search grid");

      int cellLID = locator_.findEnclosingCell(guess, x);

      int vecIndex = (*elemToVecIndexMap_)[cellLID];
      double vol = elemWeightVec_.getElement(vecIndex);
      vec.addToElement(vecIndex, particleWeight/vol);
    }
}


