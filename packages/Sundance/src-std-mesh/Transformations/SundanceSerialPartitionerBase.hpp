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

#ifndef SUNDANCE_SERIALPARTITIONERBASE_H
#define SUNDANCE_SERIALPARTITIONERBASE_H

#include "SundanceDefs.hpp"
#include "SundanceMap.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceHandle.hpp"

namespace Sundance
{
/**
 * Base class for mesh partitioners that run in serial
 */
class SerialPartitionerBase
{
public:
  
  /** */
  virtual ~SerialPartitionerBase(){;}

  /** */
  void getNeighbors(const Mesh& mesh, 
    Array<Array<int> >& neighbors, int& nEdges) const ;

  /** */
  Set<int> arrayToSet(const Array<int>& a) const ;

  /** */
  virtual void getAssignments(const Mesh& mesh, int np, 
    Array<int>& assignments) const = 0 ;

  /** */
  Array<Mesh> makeMeshParts(const Mesh& mesh, int np,
    Array<Sundance::Map<int, int> >& oldElemLIDToNewLIDMap,
    Array<Sundance::Map<int, int> >& oldVertLIDToNewLIDMap
    ) const ;

  /** */
  void getOffProcData(int p, 
    const Array<int>& elemAssignments,
    const Array<int>& nodeAssignments,
    Set<int>& offProcNodes,
    Set<int>& offProcElems) const ;

  /** 
   * 
   */
  void getNodeAssignments(int nProc, 
    const Array<int>& elemAssignments,
    Array<int>& nodeAssignments,
    Array<int>& nodeOwnerElems,
    Array<int>& nodesPerProc) const ;

  /** */
  void getElemsPerProc(int nProc, 
    const Array<int>& elemAssignments,
    Array<int>& elemsPerProc) const ;

  /** Remap global element or node 
   * numberings so that each processor owns sequentially-numbered
   * global indexes. */
  void remapEntities(const Array<int>& assignments, int nProc,
    Array<int>& entityMap) const ;


private:

  
  int max(const Set<int>& s) const ;
  mutable Array<Set<int> > elemVerts_;
  mutable Array<Set<int> > vertElems_;
};
}

#endif
