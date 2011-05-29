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
#include "SundanceUnfoldPeriodicDF.hpp"
#include "SundancePeriodicMesh1D.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundanceMeshSource.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceDOFMapBase.hpp"
#include "Teuchos_MPIContainerComm.hpp"
#include "Teuchos_TimeMonitor.hpp"



namespace Sundance
{
using namespace Teuchos;
using namespace TSFExtended;

Mesh unfoldPeriodicMesh(const Mesh& mesh)
{
  const MeshBase* mb = mesh.ptr().get();
  const PeriodicMesh1D* pm = dynamic_cast<const PeriodicMesh1D*>(mb);

  TEST_FOR_EXCEPT(pm==0);

  int numElems = mesh.numCells(1);
  double a = mesh.nodePosition(0)[0];
  double b = mesh.nodePosition(numElems)[0];

  MeshType meshType = new BasicSimplicialMeshType();
  MeshSource mesher = new PartitionedLineMesher(a, b, numElems, meshType,
    MPIComm::self());

  Mesh rtn = mesher.getMesh();

  return rtn;
}


DiscreteSpace unfoldPeriodicDiscreteSpace(const DiscreteSpace& space)
{
  VectorType<double> vecType = space.vecType();
  Mesh mesh = unfoldPeriodicMesh(space.mesh());
  BasisArray basis = space.basis();

  return DiscreteSpace(mesh, basis, vecType);
}

Expr unfoldPeriodicDiscreteFunction(const Expr& f)
{
  const DiscreteFunction* df = DiscreteFunction::discFunc(f);
  TEST_FOR_EXCEPT(df==0);


  DiscreteSpace perSpace = df->discreteSpace();
  DiscreteSpace space = unfoldPeriodicDiscreteSpace(perSpace);
  
  Mesh oldMesh = perSpace.mesh();
  Mesh newMesh = space.mesh();

  Vector<double> oldVec = df->getVector();
  Vector<double> newVec = space.createVector();

  const RCP<DOFMapBase>& oldMap = perSpace.map();
  const RCP<DOFMapBase>& newMap = space.map();

  /* Copy the element DOFs to the new vector. There are no duplicated
   * elements so this map is one-to-one */
  Array<int> oldCellLID(oldMesh.numCells(1));
  for (int i=0; i<oldMesh.numCells(1); i++)
  {
    oldCellLID[i] = i;
  }
  RCP<const Set<int> > funcs = oldMap->allowedFuncsOnCellBatch(1, oldCellLID);

  Array<Array<int> > oldElemDofs;
  Array<Array<int> > newElemDofs;
  Array<int> oldNNodes;
  RCP<const MapStructure> oldMapStruct ;
  RCP<const MapStructure> newMapStruct ;

  if (funcs->size())
  {
    oldMapStruct 
      = oldMap->getDOFsForCellBatch(1, oldCellLID, *funcs, oldElemDofs, 
        oldNNodes, 0);
    
    newMapStruct 
      = newMap->getDOFsForCellBatch(1, oldCellLID, *funcs, newElemDofs, 
        oldNNodes, 0);
    
    for (int chunk=0; chunk<oldMapStruct->numBasisChunks(); chunk++)
    {
      int nf = oldMapStruct->numFuncs(chunk);
      int nDofsPerCell = oldNNodes[chunk] * nf;
      for (int c=0; c<oldCellLID.size(); c++)
      {
        for (int d=0; d<nDofsPerCell; d++)
        {
          int oldDof = oldElemDofs[chunk][nDofsPerCell*c + d];
          int newDof = newElemDofs[chunk][nDofsPerCell*c + d];
          newVec.setElement(newDof, oldVec.getElement(oldDof));
        }
      }
    }
  }

  /* Copy the vertex dofs to the new vector. The data at v=0 are duplicated */
  Array<int> oldVertLID(oldMesh.numCells(0));
  Array<int> newVertLID(newMesh.numCells(0));
  for (int i=0; i<oldMesh.numCells(0); i++)
  {
    oldVertLID[i] = i;
  }
  for (int i=0; i<newMesh.numCells(0); i++)
  {
    newVertLID[i] = i;
  }
  if (funcs->size())
  {
    funcs = oldMap->allowedFuncsOnCellBatch(0, oldCellLID);
    
    oldMapStruct = oldMap->getDOFsForCellBatch(0, oldVertLID, *funcs, 
      oldElemDofs, 
      oldNNodes, 0);
    newMapStruct 
      = newMap->getDOFsForCellBatch(0, newVertLID, *funcs, newElemDofs, 
        oldNNodes, 0);
    
    for (int chunk=0; chunk<oldMapStruct->numBasisChunks(); chunk++)
    {
      int nf = oldMapStruct->numFuncs(chunk);
      int nDofsPerCell = oldNNodes[chunk] * nf;
      for (int c=0; c<newVertLID.size(); c++)
      {
        if (c < oldVertLID.size())
        {
          for (int d=0; d<nDofsPerCell; d++)
          {
            int oldDof = oldElemDofs[chunk][nDofsPerCell*c + d];
            int newDof = newElemDofs[chunk][nDofsPerCell*c + d];
            newVec.setElement(newDof, oldVec.getElement(oldDof));
          }
        }
        else
        {
          for (int d=0; d<nDofsPerCell; d++)
          {
            int oldDof = oldElemDofs[chunk][d];
            int newDof = newElemDofs[chunk][nDofsPerCell*c + d];
            newVec.setElement(newDof, oldVec.getElement(oldDof));
          }
        }
      }
    }
  }

  return new DiscreteFunction(space, newVec);
}


}


