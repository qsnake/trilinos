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

#ifndef SUNDANCE_INCREMENTALLYCREATABLEMESH_H
#define SUNDANCE_INCREMENTALLYCREATABLEMESH_H



#include "SundanceDefs.hpp"
#include "SundanceMeshBase.hpp"
#include "Teuchos_Array.hpp"

namespace Sundance
{
/**
 * IncrementallyCreatableMesh is an interface for the creation of meshes
 * by adding vertices, then elements.
 */
class IncrementallyCreatableMesh : public MeshBase
{
public:
  /** Construct an empty mesh of the given dimension, distributed over
   * processors in the MPI communicator*/
  IncrementallyCreatableMesh(int dim, const MPIComm& comm, 
    const MeshEntityOrder& meshOrder)
    : MeshBase(dim, comm, meshOrder) {;}

  /** virtual dtor */
  virtual ~IncrementallyCreatableMesh(){;}

  /** Optional preallocation of space for an estimated number of vertices */
  virtual void estimateNumVertices(int nPts) {;}

  /** Optional preallocation of space for an estimated number of elements */
  virtual void estimateNumElements(int nElems) {;}

      

  /** 
   * Add new new vertex to the mesh.
   * \param globalIndex the GID of the new vertex
   * \param x the spatial position of the new vertex
   * \param ownerProcID the processor that "owns" this vertex 
   * \param label a label for this vertex (optionally used in setting loads, boundary
   * conditions, etc)
   * \return the LID of the vertex.
   */
  virtual int addVertex(int globalIndex, const Point& x,
    int ownerProcID, int label) = 0 ;

  /** 
   * Add a new element to the mesh.
   * \param globalIndex the GID of the new element
   * \param vertexGIDs tuple of GIDs for the vertices defining this element. 
   * \param ownerProcID the processor that "owns" this element
   * \param label a label for this element (optionally used in setting loads, 
   * material properties, etc)
   * \return the LID of the element
   */
  virtual int addElement(int globalIndex, const Array<int>& vertexGIDs,
    int ownerProcID, int label) = 0 ;
  

};

}



#endif
