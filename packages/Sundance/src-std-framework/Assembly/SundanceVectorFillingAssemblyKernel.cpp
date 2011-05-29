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
#include "SundanceVectorFillingAssemblyKernel.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "TSFLoadableBlockVector.hpp"



using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using namespace TSFExtended;
using std::setw;
using std::endl;
      
static Time& vecInsertTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("vector insertion"); 
  return *rtn;
}

VectorFillingAssemblyKernel::VectorFillingAssemblyKernel(
  const Array<RCP<DOFMapBase> >& dofMap,
  const Array<RCP<Array<int> > >& isBCIndex,
  const Array<int>& lowestLocalIndex,
  Array<Vector<double> >& b,
  bool partitionBCs,
  int verbosity
  )
  : AssemblyKernelBase(verbosity),
    b_(b),
    vec_(b.size()),
    mapBundle_(dofMap, isBCIndex, lowestLocalIndex, partitionBCs, verbosity)
{
  Tabs tab0;

  SUNDANCE_MSG1(verb(), tab0 << "VectorFillingAssemblyKernel ctor");
  
  int numBlocks = dofMap.size();

  for (int i=0; i<b_.size(); i++)
  {
    vec_[i].resize(numBlocks);
    for (int block=0; block<numBlocks; block++)
    {
      Tabs tab1;
      SUNDANCE_MSG1(verb(), tab1 << "getting vector for block b=" 
        << block << " of " << numBlocks);
      Vector<double> vecBlock; 
      if (partitionBCs && numBlocks == 1)
      {
        Tabs tab2;
        SUNDANCE_MSG1(verb(), tab2 << "making loadable block vector");
        vecBlock = b[i];
        int lowestRow = mapBundle_.lowestLocalIndex(block) ;
        int highestRow = lowestRow + mapBundle_.dofMap(block)->numLocalDOFs();
        vec_[i][block] = 
          rcp(new LoadableBlockVector(vecBlock, lowestRow,
              highestRow, mapBundle_.isBCIndex(block)));
      }
      else
      {
        vecBlock = b[i].getBlock(block);
        vec_[i][block] = rcp_dynamic_cast<LoadableVector<double> >(vecBlock.ptr());
      }
      TEST_FOR_EXCEPTION(vec_[i][block].get()==0, RuntimeError,
        "vector block " << block << " is not loadable");
      vecBlock.zero();
    }
  }
  SUNDANCE_MSG1(verb(), tab0 << "done VectorFillingAssemblyKernel ctor");
}


void VectorFillingAssemblyKernel::buildLocalDOFMaps(
  const RCP<StdFwkEvalMediator>& mediator,
  IntegrationCellSpecifier intCellSpec,
  const Array<Set<int> >& requiredFuncs) 
{
  mapBundle_.buildLocalDOFMaps(mediator, intCellSpec, requiredFuncs,
    verb());
}


void VectorFillingAssemblyKernel::insertLocalVectorBatch(
  bool isBCRqc,
  bool useCofacetCells,
  const Array<int>& funcID,  
  const Array<int>& funcBlock, 
  const Array<int>& mvIndices, 
  const Array<double>& localValues) const
{
  TimeMonitor timer(vecInsertTimer());
  Tabs tab0;

  SUNDANCE_MSG1(verb(), tab0 << "inserting local vector batch");
  SUNDANCE_MSG4(verb(), tab0 << "vector values are " << localValues);

  const MapBundle& mb = mapBundle_;
  int nCells = mb.nCells();

  for (int i=0; i<funcID.size(); i++)
  {
    Tabs tab1;
    SUNDANCE_MSG2(verb(), tab1 << "function ID = "<< funcID[i] 
      << " of " << funcID.size());
    SUNDANCE_MSG2(verb(), tab1 << "is BC eqn = " << isBCRqc);
    SUNDANCE_MSG2(verb(), tab1 << "num cells = " << nCells);
    SUNDANCE_MSG2(verb(), tab1 << "using cofacet cells = " << useCofacetCells);
    SUNDANCE_MSG2(verb(), tab1 << "multivector index = " 
      << mvIndices[i]);

    /* First, find the block associated with the current function
     * so that we can find the appropriate DOF information */
    int block = funcBlock[i];

    const RCP<DOFMapBase>& dofMap = mb.dofMap(block);
    int lowestLocalRow = mb.lowestLocalIndex(block);

    int chunk = mb.mapStruct(block, useCofacetCells)->chunkForFuncID(funcID[i]);
    SUNDANCE_MSG2(verb(), tab1 << "chunk = " << chunk);

    int funcIndex = mb.mapStruct(block, useCofacetCells)->indexForFuncID(funcID[i]);
    SUNDANCE_MSG2(verb(), tab1 << "func offset into local DOF map = " 
      << funcIndex);

    const Array<int>& dofs = mb.localDOFs(block, useCofacetCells, chunk);
    SUNDANCE_MSG4(verb(), tab1 << "local dofs = " << dofs);    

    int nFuncs = mb.mapStruct(block, useCofacetCells)->numFuncs(chunk);
    SUNDANCE_MSG2(verb(), tab1 << "num funcs in chunk = " << nFuncs);

    int nNodes = mb.nNodesInChunk(block, useCofacetCells, chunk);
    SUNDANCE_MSG2(verb(), tab1 << "num nodes in chunk = " << nNodes);

    const Array<int>& isBCIndex = *(mb.isBCIndex(block));

    /* At this point, we can start to load the elements */
    int r=0;
    RCP<TSFExtended::LoadableVector<double> > vecBlock 
      = vec_[mvIndices[i]][block];

    FancyOStream& os = Out::os();

    for (int c=0; c<nCells; c++)
    {
      Tabs tab2;
      SUNDANCE_MSG2(verb(), tab2 << "cell = " << c << " of " << nCells);
      for (int n=0; n<nNodes; n++, r++)
      {
        Tabs tab3;
        int rowIndex = dofs[(c*nFuncs + funcIndex)*nNodes + n];
        int localRowIndex = rowIndex - lowestLocalRow;
        if (verb() >= 2) os << tab3 << "n=" << setw(4) << n 
                            << " G=" << setw(8) << rowIndex 
                            << " L=" << setw(8) << localRowIndex;
        if (!(dofMap->isLocalDOF(rowIndex))
          || isBCRqc!=isBCIndex[localRowIndex]) 
        {
          if (verb() >= 2)
          {
            if (!dofMap->isLocalDOF(rowIndex)) 
            {
              os << " --- skipping (is non-local) ---" << std::endl;
            }
            else if (!isBCRqc && isBCIndex[localRowIndex])
            {
              os << " --- skipping (is BC row) ---" << std::endl;
            }
            else
            {
              os << " --- skipping (is non-BC row) ---" << std::endl;
            }
          }
        }
        else
        {
          if (verb() >= 2) os << setw(15) << localValues[r] << std::endl;
          vecBlock->addToElement(rowIndex, localValues[r]);
        }
      }
    }
  }
  SUNDANCE_MSG1(verb(), tab0 << "...done vector insertion");
}


