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

#include "SundanceStdFwkEvalMediator.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceExceptions.hpp"


using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;





StdFwkEvalMediator::StdFwkEvalMediator(const Mesh& mesh, int cellDim)
  : AbstractEvalMediator(),
    mesh_(mesh),
    cellDim_(cellDim),
    cellType_(NullCell),
    maxCellType_(NullCell),
    isInternalBdry_(false),
    forbidCofacetIntegrations_(false),
    cellLID_(),
    JVol_(rcp(new CellJacobianBatch())),
    JTrans_(rcp(new CellJacobianBatch())),
    facetIndices_(rcp(new Array<int>())),
    maxCellLIDs_(rcp(new Array<int>())),
    cofacetCellsAreReady_(false),
    cacheIsValid_(false),
    jCacheIsValid_(false),
    fCache_(),
    dfCache_(),
    localValueCache_(),
    mapStructCache_(),
    fCacheIsValid_(),
    dfCacheIsValid_(),
    localValueCacheIsValid_()
{;}

void StdFwkEvalMediator::setCellType(const CellType& cellType,
  const CellType& maxCellType,
  bool isInternalBdry) 
{
  cellType_=cellType; 
  cellDim_ = dimension(cellType);
  maxCellType_ = maxCellType;
  isInternalBdry_ = isInternalBdry;
  cacheIsValid() = false; 
  jCacheIsValid_=false;
  cofacetCellsAreReady_ = false;
  forbidCofacetIntegrations_ = isInternalBdry_ ;
//    && !ElementIntegral::alwaysUseCofacets();
}


void StdFwkEvalMediator::setIntegrationSpec(
  IntegrationCellSpecifier intCellSpec)
{
  intCellSpec_ = intCellSpec;
}

void StdFwkEvalMediator::setCellBatch(
  const RCP<const Array<int> >& cellLID) 
{
  Tabs tab(0);
  SUNDANCE_MSG1(verb(), tab << "eval med setting cell batch");
  Tabs tab1;

  cellLID_ = cellLID; 
  cacheIsValid() = false; 
  jCacheIsValid_=false;
  cofacetCellsAreReady_ = false;

  SUNDANCE_MSG2(verb(), tab1 << "getting volume Jacobians");
  mesh_.getJacobians(cellDim(), *cellLID, *JVol_);
  if (intCellSpec_!=NoTermsNeedCofacets) setupFacetTransformations();

  /* mark the function caches as invalid */
  SUNDANCE_MSG2(verb(), tab1 << "flushing old function caches");
  Map<const DiscreteFunctionData*, bool>::iterator iter;
  for (iter = fCacheIsValid_.begin(); iter != fCacheIsValid_.end(); iter++)
    {
      iter->second = false;
    }
  for (iter = dfCacheIsValid_.begin(); iter != dfCacheIsValid_.end(); iter++)
    {
      iter->second = false;
    }
  for (iter = localValueCacheIsValid_.begin(); iter != localValueCacheIsValid_.end(); iter++)
    {
      iter->second = false;
    }
}

void StdFwkEvalMediator::setupFacetTransformations() const 
{
  Tabs tab(0);
  SUNDANCE_MSG1(verb(), tab << "setting up facet transformations");
  Tabs tab1;

  const Array<int>& cells = *cellLID_;

  SUNDANCE_MSG2(verb(), tab1 << "num cells in batch = " << cells.size());
  SUNDANCE_MSG2(verb(), tab1 << "cell dim = " << cellDim());
  SUNDANCE_MSG2(verb(), tab1 << "num d-cells in mesh = " << mesh_.numCells(cellDim()));
  

  facetIndices_->resize(cells.size());
  maxCellLIDs_->resize(cells.size());
  cofacetCellsAreReady_ = true;

  for (int c=0; c<cells.size(); c++)
    {
      Tabs tab2;
      SUNDANCE_MSG2(verb(), tab2 << "c = " << c);
      SUNDANCE_MSG2(verb(), tab2 << "cells[c] = " << cells[c]);
      (*maxCellLIDs_)[c] 
        = mesh_.maxCofacetLID(cellDim(), cells[c], 0, (*facetIndices_)[c]);
    }

  SUNDANCE_MSG2(verb(), tab1 << "getting facet Jacobians");
  mesh_.getJacobians(mesh_.spatialDim(), *maxCellLIDs_, *JTrans_);
  SUNDANCE_MSG2(verb(), tab << "setting up facet transformations");
}



const CellJacobianBatch& StdFwkEvalMediator::JTrans() const
{
  /* If we're integrating a derivative on a boundary, JVol and JTrans will be
   * different. Otherwise, they'll be the same, and we use JVol for both
   * volume computations and vector transformations */
  if (intCellSpec_ != NoTermsNeedCofacets) 
  {
    return *JTrans_;
  }
  else
  {
    return *JVol_;
  }
}
