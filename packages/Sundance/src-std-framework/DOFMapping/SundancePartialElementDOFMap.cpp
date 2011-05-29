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

#include "SundanceMap.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"
#include "SundanceOrderedTuple.hpp"
#include "SundancePartialElementDOFMap.hpp"
#include "SundanceLagrange.hpp"
#include "Teuchos_MPIContainerComm.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace Sundance;
using namespace Teuchos;

PartialElementDOFMap::PartialElementDOFMap(const Mesh& mesh, 
  const CellFilter& subdomain,
  int nFuncs,
  int setupVerb)
  : DOFMapBase(mesh, setupVerb),
    dim_(mesh.spatialDim()),
    nFuncs_(nFuncs),
    nElems_(mesh.numCells(mesh.spatialDim())),
    subdomain_(subdomain),
    funcDomains_(nFuncs, subdomain),
    elemDofs_(nElems_ * nFuncs, -1),
    structure_(rcp(new MapStructure(nFuncs_, rcp(new Lagrange(0))))),
    allFuncs_()
{
  init();
  Set<int> tmp;
  for (int i=0; i<nFuncs_; i++) tmp.put(i);
  allFuncs_ = rcp(new Set<int>(tmp));
}


RCP<const MapStructure> 
PartialElementDOFMap::getDOFsForCellBatch(int cellDim,
  const Array<int>& cellLID,
  const Set<int>& requestedFuncSet,
  Array<Array<int> >& dofs,
  Array<int>& nNodes,
  int verbosity) const
{
  TimeMonitor timer(batchedDofLookupTimer());


  Tabs tab;
  SUNDANCE_MSG3(verbosity, 
    tab << "PartialElementDOFMap::getDOFsForCellBatch(): cellDim=" 
    << cellDim
    << " cellLID=" << cellLID);


  dofs.resize(1);
  nNodes.resize(1);
  nNodes[0] = 1;

  int nCells = cellLID.size();


  if (cellDim == dim_)
  {
    dofs[0].resize(nCells * nFuncs_);
    Array<int>& dof0 = dofs[0];
      
    for (int c=0; c<nCells; c++)
    {
      for (int i=0; i<nFuncs_; i++)
      {
        dof0[c*nFuncs_ + i] = elemDofs_[cellLID[c]*nFuncs_+i];
      }
    }
  }
  else
  {
    dofs[0].resize(nCells * nFuncs_);
    Array<int> cofacetLIDs(nCells);
      
    for (int c=0; c<nCells; c++)
    {
      TEST_FOR_EXCEPTION(mesh().numMaxCofacets(cellDim, cellLID[c]) > 1,
        RuntimeError,
        "Attempt to do a trace of a L0 basis on a "
        "lower-dimensional cell having more than one "
        "maximal cofacets");
      int myFacetIndex = -1;
      cofacetLIDs[c] = mesh().maxCofacetLID(cellDim, cellLID[c],
        0, myFacetIndex);
    }

    getDOFsForCellBatch(dim_, cofacetLIDs, requestedFuncSet, dofs,
      nNodes, verbosity);
  }

  return structure_;
}



void PartialElementDOFMap::init() 
{ 
  Tabs tab;

  SUNDANCE_MSG1(setupVerb(), tab << "initializing element DOF map");

  int cellDim = mesh().spatialDim();

  Array<Array<int> > remoteElems(mesh().comm().getNProc());

  /* start the DOF count at zero. */
  int nextDOF = 0;
  int owner;

  /* Assign node DOFs for locally owned maximal cells */
  CellSet maxCells = subdomain_.getCells(mesh());

  for (CellIterator iter=maxCells.begin(); iter != maxCells.end(); iter++)
  {
    int cellLID = *iter;
    /* the cell may be owned by another processor */
    if (isRemote(cellDim, cellLID, owner))
    {
      int elemGID = mesh().mapLIDToGID(cellDim, cellLID);
      remoteElems[owner].append(elemGID);
    }                  
    else /* we can assign a DOF locally */
    {
      /* assign DOFs */
      for (int i=0; i<nFuncs_; i++)
      {
        elemDofs_[cellLID*nFuncs_ + i] = nextDOF;
        nextDOF++;
      }
    }
  }
  
  /* Compute offsets for each processor */
  int localCount = nextDOF;
  computeOffsets(localCount);
  
  /* Resolve remote DOF numbers */
  shareRemoteDOFs(remoteElems);
}


RCP<const Set<int> > PartialElementDOFMap
::allowedFuncsOnCellBatch(int cellDim,
  const Array<int>& cellLID) const 
{
  static RCP<const Set<int> > empty = rcp(new Set<int>());

  if (cellDim != dim_) return empty;
  for (int i=0; i<cellLID.size(); i++)
  {
    if (elemDofs_[cellLID[i]*nFuncs_] < 0) return empty;
  }
  return allFuncs_;
}


void PartialElementDOFMap::computeOffsets(int localCount)
{
  Array<int> dofOffsets;
  int totalDOFCount;
  int myOffset = 0;
  int np = mesh().comm().getNProc();
  if (np > 1)
  {
    MPIContainerComm<int>::accumulate(localCount, dofOffsets, totalDOFCount,
      mesh().comm());
    myOffset = dofOffsets[mesh().comm().getRank()];

    int nDofs = nElems_ * nFuncs_;
    for (int i=0; i<nDofs; i++)
    {
      if (elemDofs_[i] >= 0) elemDofs_[i] += myOffset;
    }
  }
  else
  {
    totalDOFCount = localCount;
  }
  
  setLowestLocalDOF(myOffset);
  setNumLocalDOFs(localCount);
  setTotalNumDOFs(totalDOFCount);
}


void PartialElementDOFMap::shareRemoteDOFs(const Array<Array<int> >& outgoingCellRequests)
{
  int cellDim = mesh().spatialDim();

  int np = mesh().comm().getNProc();
  if (np==1) return;
  int rank = mesh().comm().getRank();

  Array<Array<int> > incomingCellRequests;
  Array<Array<int> > outgoingDOFs(np);
  Array<Array<int> > incomingDOFs;

  SUNDANCE_MSG2(setupVerb(),  
    "p=" << mesh().comm().getRank()
    << "synchronizing DOFs for cells of dimension 0");
  SUNDANCE_MSG2(setupVerb(),  
    "p=" << mesh().comm().getRank()
    << " sending cell reqs d=0, GID=" 
    << outgoingCellRequests);

  /* share the cell requests */
  MPIContainerComm<int>::allToAll(outgoingCellRequests, 
    incomingCellRequests,
    mesh().comm());
  
  /* get DOF numbers for the zeroth function index on every node that's been 
   * requested by someone else */
  for (int p=0; p<np; p++)
  {
    if (p==rank) continue;
    const Array<int>& requestsFromProc = incomingCellRequests[p];
    int nReq = requestsFromProc.size();

    SUNDANCE_MSG4(setupVerb(), "p=" << mesh().comm().getRank() 
      << " recv'd from proc=" << p
      << " reqs for DOFs for cells " 
      << requestsFromProc);

    outgoingDOFs[p].resize(nReq);

    for (int c=0; c<nReq; c++)
    {
      int GID = requestsFromProc[c];
      SUNDANCE_MSG3(setupVerb(),  
        "p=" << rank
        << " processing zero-cell with GID=" << GID); 
      int LID = mesh().mapGIDToLID(cellDim, GID);
      SUNDANCE_MSG3(setupVerb(),  
        "p=" << rank
        << " LID=" << LID << " dofs=" 
        << elemDofs_[LID*nFuncs_]);
      outgoingDOFs[p][c] = elemDofs_[LID*nFuncs_];
      SUNDANCE_MSG3(setupVerb(),  
        "p=" << rank
        << " done processing cell with GID=" << GID);
    }
  }

  SUNDANCE_MSG2(setupVerb(),  
    "p=" << mesh().comm().getRank()
    << "answering DOF requests for cells of dimension 0");

  /* share the DOF numbers */
  MPIContainerComm<int>::allToAll(outgoingDOFs,
    incomingDOFs,
    mesh().comm());

  SUNDANCE_MSG2(setupVerb(),  
    "p=" << mesh().comm().getRank()
    << "communicated DOF answers for cells of dimension 0" );

  
  /* now assign the DOFs from the other procs */

  for (int p=0; p<mesh().comm().getNProc(); p++)
  {
    if (p==mesh().comm().getRank()) continue;
    const Array<int>& dofsFromProc = incomingDOFs[p];
    int numCells = dofsFromProc.size();
    for (int c=0; c<numCells; c++)
    {
      int cellGID = outgoingCellRequests[p][c];
      int cellLID = mesh().mapGIDToLID(cellDim, cellGID);
      int dof = dofsFromProc[c];
      for (int i=0; i<nFuncs_; i++)
      {
        elemDofs_[cellLID*nFuncs_ + i] = dof+i;
        addGhostIndex(dof+i);
      }
    }
  }
}


void PartialElementDOFMap::print(std::ostream& os) const
{
  os << elemDofs_ << std::endl;
}
