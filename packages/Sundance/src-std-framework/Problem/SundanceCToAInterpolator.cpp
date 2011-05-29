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

#include "SundanceCToAInterpolator.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceLagrange.hpp"
#include "SundanceDiscreteFuncElement.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using namespace TSFExtended;


static Time& particleInterpolationTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("particle interpolation"); 
  return *rtn;
}

CToAInterpolator::CToAInterpolator(const AToCPointLocator& locator,
  const Expr& field)
  : dim_(locator.mesh().spatialDim()),
    nFacets_(dim_+1),
    rangeDim_(-1),
    elemToVecValuesMap_(),
    locator_(locator)
{
  updateField(field);
}

void CToAInterpolator::updateField(const Expr& field)
{
  int newRangeDim = field.size();
  if (newRangeDim != rangeDim_)
  {
    rangeDim_ = newRangeDim;
    elemToVecValuesMap_ 
      = rcp(new Array<double>(rangeDim_ * locator_.mesh().numCells(dim_) * nFacets_));
  }

  int nCells = locator_.mesh().numCells(dim_);

  const DiscreteFunction* df = DiscreteFunction::discFunc(field);
  TEST_FOR_EXCEPTION(df == 0, RuntimeError,
    "discrete function expected in "
    "CToAInterpolator::updateField()");
  
  Vector<double> vec = df->getVector();      

  Array<int> cellLID(nCells);
  for (int i=0; i<nCells; i++)
  {
    cellLID[i] = i;
  }
  Array<Array<int> > dofs;
  Array<int> nNodes;
  Set<int> funcs;
  for (int i=0; i<rangeDim_; i++) funcs.put(i);

  
  df->map()->getDOFsForCellBatch(dim_, cellLID, funcs, dofs, nNodes,0);

  const Array<int>& dofs0 = dofs[0];

  for (int c=0; c<cellLID.size(); c++)
  {
    for (int n=0; n<nFacets_; n++)
    {
      for (int f=0; f<rangeDim_; f++)
      {
        int dof = dofs0[(c*rangeDim_ + f)*nFacets_ + n];
        (*elemToVecValuesMap_)[(cellLID[c]*nFacets_ + n)*rangeDim_ + f]
          = vec.getElement(dof);
      }
    }
  }
}


void CToAInterpolator::interpolate(const Teuchos::Array<double>& positions,
  Teuchos::Array<double>& results) const
{
  TimeMonitor timer(particleInterpolationTimer());

  TEST_FOR_EXCEPTION(positions.size() % dim_ != 0, RuntimeError,
    "vector of coordinates should by an integer multiple "
    "of the spatial dimension");

  int nPts = positions.size() / dim_;

  results.resize(rangeDim_ * nPts);

  Array<double> xLocal(dim_);

  for (int i=0; i<nPts; i++)
  {
    const double* x = &(positions[dim_*i]);

    int guess = locator_.guessCell(x);

    TEST_FOR_EXCEPTION(guess < 0, RuntimeError, "particle position "
      << x << " out of search grid");

      
    int cellLID = locator_.findEnclosingCell(guess, x, &(xLocal[0]));

    if (dim_==2)
    {
      double s = xLocal[0];
      double t = xLocal[1];
      Array<double> phi(nFacets_);
      phi[0] = 1.0-s-t;
      phi[1] = s;
      phi[2] = t;
      for (int f=0; f<rangeDim_; f++) results[rangeDim_*i+f] = 0.0;
      for (int n=0; n<nFacets_; n++)
      {
        for (int f=0; f<rangeDim_; f++)
        {
          results[rangeDim_*i+f] += phi[n]*(*elemToVecValuesMap_)[(cellLID*nFacets_ + n)*rangeDim_ + f];
        }
      }
    }
  }

}

