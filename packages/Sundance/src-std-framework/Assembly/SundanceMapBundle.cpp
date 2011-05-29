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

#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceMapBundle.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceStdFwkEvalMediator.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

      

MapBundle::MapBundle(
  const Array<RCP<DOFMapBase> >& dofMap,
  const Array<RCP<Array<int> > >& isBCIndex,
  const Array<int>& lowestLocalIndex,
  bool partitionBCs,
  int verb
  )
  : verb_(verb),
    dofMap_(dofMap),
    isBCIndex_(isBCIndex),
    lowestLocalIndex_(lowestLocalIndex),
    localDOFMap_(rcp(new LocalDOFMap(dofMap_.size(), verb))),
    cofacetLocalDOFMap_(rcp(new LocalDOFMap(dofMap_.size(), verb)))
{}

int MapBundle::nCells() const 
{
  TEST_FOR_EXCEPTION(localDOFMap_->isUnused() 
    && cofacetLocalDOFMap_->isUnused(), RuntimeError,
    "no local DOF maps defined in call to MapBundle::nCells()");

  
  if (cofacetLocalDOFMap_->isUnused())
  {
    return localDOFMap_->nCells();
  }
  else if (localDOFMap_->isUnused())
  {
    return cofacetLocalDOFMap_->nCells();
  }
  else
  {
    TEST_FOR_EXCEPTION(localDOFMap_->nCells() != cofacetLocalDOFMap_->nCells(),
      RuntimeError,
      "mismatched cell counts in MapBundle::nCells()");
    return cofacetLocalDOFMap_->nCells();
  }
}

RCP<const Array<int> > MapBundle::workSet(int block,
  bool useCofacets) const
{
  return chooseMap(block, useCofacets)->cellLIDs();
}


const RCP<LocalDOFMap>& MapBundle::chooseMap(
  int block, bool useCofacets) const
{
  if (useCofacets)
  {
    TEST_FOR_EXCEPTION(cofacetLocalDOFMap_->isUnused(block),
      RuntimeError,
      "request for unavailable cofacet-based local map for block = " << block);
    return cofacetLocalDOFMap_;
  }
  else
  {
    TEST_FOR_EXCEPTION(localDOFMap_->isUnused(block),
      RuntimeError,
      "request for unavailable local map for block = " << block);
    return localDOFMap_;
  }
}




void MapBundle::buildLocalDOFMaps(
  const RCP<StdFwkEvalMediator>& mediator,
  IntegrationCellSpecifier intCellSpec,
  const Array<Set<int> >& requiredFuncs,
  int verbosity)
{
  Tabs tab;

  int numBlocks = dofMap_.size();

  localDOFMap_->markAsUnused();
  cofacetLocalDOFMap_->markAsUnused();
  localDOFMap_->setVerbosity(verbosity);
  cofacetLocalDOFMap_->setVerbosity(verbosity);

  int maxCellDim = mediator->maxCellDim();
  int cellDim = mediator->cellDim();

  SUNDANCE_MSG3(verbosity, tab << "cell dim=" << cellDim);
  SUNDANCE_MSG3(verbosity, tab << "max cell dim=" << maxCellDim);

  for (int b=0; b<numBlocks; b++)
  {   
    Tabs tab2;
    SUNDANCE_MSG3(verbosity, tab2 << "getting dofs for block " 
      << b << " of " << numBlocks);
        
    if (intCellSpec != AllTermsNeedCofacets)
    {
      Tabs tab3;
      SUNDANCE_MSG3(verbosity, tab3 << "getting ordinary dofs");

      if (!localDOFMap_->hasCells()) 
      {
        SUNDANCE_MSG3(verbosity, tab3 << "setting cells of dim " 
          << cellDim);
        localDOFMap_->setCells(cellDim, maxCellDim, mediator->cellLID());
      }     
      localDOFMap_->fillBlock(b, dofMap_[b], requiredFuncs);
      localDOFMap_->markAsUsed(b);
    }
    else
    {
      Tabs tab3;
      SUNDANCE_MSG3(verbosity, tab3 << "ordinary dofs not needed for block " << b);
    }

        
    if (intCellSpec != NoTermsNeedCofacets)
    {
      Tabs tab3;
      SUNDANCE_MSG3(verbosity, tab3 << "getting cofacet dofs");
      SUNDANCE_MSG3(verbosity, tab3 << "cofacet cells="
        << *(mediator->cofacetCellLID()));

      if (!cofacetLocalDOFMap_->hasCells()) 
        cofacetLocalDOFMap_->setCells(maxCellDim, 
          maxCellDim, mediator->cofacetCellLID());

      cofacetLocalDOFMap_->fillBlock(b, dofMap_[b], requiredFuncs);
      cofacetLocalDOFMap_->markAsUsed(b);
    }
    else
    {
      Tabs tab3;
      SUNDANCE_MSG3(verbosity, tab3 << "cofacet dofs not needed for block " << b);
    }
        
  }

  if (intCellSpec != AllTermsNeedCofacets)
  {
    SUNDANCE_MSG4(verbosity, tab << "local DOF values " << *localDOFMap_);
  }

  if (intCellSpec != NoTermsNeedCofacets)
  {
    SUNDANCE_MSG4(verbosity, tab << "local cofacet DOF values " 
      << *cofacetLocalDOFMap_);
  }
}
