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

#ifndef SUNDANCE_MESHTYPE_H
#define SUNDANCE_MESHTYPE_H

#include "SundanceDefs.hpp"
#include "SundanceMeshTypeBase.hpp"
#include "SundanceMesh.hpp"
#include "SundanceHandle.hpp"

namespace Sundance
{
/**
 * Class MeshType is a user-level object for specification of which
 * internal mesh representation is to be used when building or reading
 * a mesh. An example of using a MeshType to control the creation 
 * of a mesh with a TriangleMeshReader is as follows: 
 * \code
 * MeshType meshType = new BasicSimplicialMeshType();
 * MeshSource meshSrc = new TriangleMeshReader("meshFile", meshType, MPIComm::world());
 * \endcode
 * The internal representation of the mesh will be as a BasicSimplicialMesh
 * object. 
 */
class MeshType : public Sundance::Handle<MeshTypeBase>
{
public:
  /** Construct an empty mesh type object */
  MeshType();

  /** Construct from a raw pointer to a mesh type subtype */
  MeshType(Sundance::Handleable<MeshTypeBase>* rawPtr);

  /** Construct from a smart pointer to a mesh type subtype */
  MeshType(const RCP<MeshTypeBase>& smartPtr);

  /** Create a mesh of the given dimension */
  Mesh createEmptyMesh(int dim, const MPIComm& comm) const ;
    
};
}

#endif
