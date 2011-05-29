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
#include "SundanceDOFMapBase.hpp"
#include "SundanceLocalDOFMap.hpp"


using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;


using std::endl;
using std::setw;


LocalDOFMap::LocalDOFMap(int numBlocks, int verb)
  : verb_(verb),
    isUsed_(numBlocks),
    hasCells_(false),
    nLocalNodesPerChunk_(rcp(new Array<Array<int> >(numBlocks))),
    mapStruct_(rcp(new Array<RCP<const MapStructure> >(numBlocks))),
    localDOFs_(rcp(new Array<Array<Array<int> > >(numBlocks))),
    cellLID_(),
    activeCellDim_(-1),
    maxCellDim_(-1)
{}

int LocalDOFMap::nCells() const 
{
  TEST_FOR_EXCEPTION(!hasCells(), RuntimeError,
    "cells not valid when LocalDOFMap::nCells() called");
  return cellLID_->size();
}

void  LocalDOFMap::markAsUnused() 
{
  for (int b=0; b<numBlocks(); b++) 
  {
    isUsed_[b] = 0;
  }
  hasCells_ = false;
}


bool  LocalDOFMap::isUnused() const 
{
  for (int b=0; b<numBlocks(); b++) 
  {
    if (isUsed_[b]) return false;
  }
  return true;
}

void LocalDOFMap::verifyValidBlock(int b) const 
{
  TEST_FOR_EXCEPTION(b < 0 || b>=numBlocks(), RuntimeError,
    "block index " << b << " out of range [0," << numBlocks() << ")");
}

void LocalDOFMap::setCells(
  int activeCellDim, 
  int maxCellDim,
  const RCP<const Array<int> >& cellLID)
{
  activeCellDim_ = activeCellDim;
  maxCellDim_ = maxCellDim_;
  cellLID_ = cellLID;
  hasCells_=true;
}


void LocalDOFMap::fillBlock(int b, const RCP<DOFMapBase>& globalMap,
  const Array<Set<int> >& requiredFuncs)
{
  Tabs tab;
  verifyValidBlock(b);

  SUNDANCE_MSG2(verb_, tab << "getting local DOFs for block=" << b 
    << ", active cell dim="
    << activeCellDim_);

  mapStruct(b) = globalMap->getDOFsForCellBatch(
    activeCellDim_,
    *cellLID_,
    requiredFuncs[b],
    localDOFs(b),
    nLocalNodesPerChunk(b),
    verb_);
}

std::ostream& LocalDOFMap::print(std::ostream& os) const
{
  Tabs tab0;
  bool noneAreUsed=true;
  for (int b=0;  b<numBlocks(); b++) 
  {
    if (isUsed(b)) {noneAreUsed = false; break;}
  }
  if (noneAreUsed)
  {
    os << std::endl << tab0 << "[empty local DOF map]";
  }
  else
  {
    os << std::endl << tab0 << "LocalDOFMap[" << std::endl;
    Tabs tab;
    os << tab << "num cells=" << cellLID_->size() << ", active cellDim="
       << activeCellDim_ << ", maxCellDim=" << maxCellDim_ << std::endl;
    os << tab << "num blocks=" << numBlocks() << std::endl << std::endl;

    for (int b=0; b<numBlocks(); b++)
    {
      Tabs tab1;
      os << tab1 << "block " << b << " of " << numBlocks()
         << *mapStruct(b) << std::endl;
      const Array<Array<int> >& dofs = localDOFs(b);
      int nChunks = mapStruct(b)->numBasisChunks();
      
      for (int c=0; c<cellLID_->size(); c++)
      {
        Tabs tab2;
        os << tab2 << "cell LID=" << (*cellLID_)[c] << std::endl;
        for (int chunk=0; chunk<nChunks; chunk++)
        {
          Tabs tab3;
          int nFuncs = mapStruct(b)->numFuncs(chunk);
          const Array<int>& funcs = mapStruct(b)->funcs(chunk);
          int nNodes = nLocalNodesPerChunk(b)[chunk];
          for (int f=0; f<nFuncs; f++)
          {
            os << tab3 << "fid=" << funcs[f] << " dofs=";
            int funcOffset = mapStruct(b)->indexForFuncID(funcs[f]);
            for (int n=0; n<nNodes; n++)
            {
              int dof = dofs[chunk][(c*nFuncs+funcOffset)*nNodes+n];
              os << setw(6) << dof;
              if (n < nNodes-1) os << ", ";
            }
            os << std::endl;
          }
        }
      }
      
    }
    os << tab0 << "] ### End LocalDOFMap" << std::endl;
  }
  return os;
}
