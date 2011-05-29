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

#ifndef SUNDANCE_MESHFILTER_H
#define SUNDANCE_MESHFILTER_H

#include "SundanceDefs.hpp"
#include "SundanceMeshTransformationBase.hpp"
#include "SundanceHandle.hpp"

namespace Sundance
{
  /**
   * MeshTransformation is the user-level interface for mesh filters, i.e.,
   * objects that take an input mesh and produce a new mesh. Examples
   * of filter operations are refinement, load balancing,
   * and extrusion from 2D to 3D. 
   *
   * <h4> Example: </h4> extrude a 2D mesh into 2D
   * \code
   * // create a 2D mesh 
   * MeshType meshType = new BasicSimplicialMeshType();
   * MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, 10, 1,
   *                                                    0.0, 1.0, 10, 1,
   *                                                    meshType);
   * Mesh mesh2D = mesher.getMesh();
   * // create a filter for extruding 2 levels between z=0.0 and z=0.2
   * MeshTransformation extruder = new ExtrusionMeshTransformation(0.0, 0.2, 2);
   * // perform the extrusion
   * Mesh mesh3D = extruder.apply(mesh2D);
   * \endcode
   */
  class MeshTransformation : public Sundance::Handle<MeshTransformationBase>
  {
  public:
    /** Construct an empty mesh filter object */
    MeshTransformation();

    /** Construct from a raw pointer to a mesh filter subtype */
    MeshTransformation(Sundance::Handleable<MeshTransformationBase>* rawPtr);

    /** Construct from a smart pointer to a mesh filter subtype */
    MeshTransformation(const RCP<MeshTransformationBase>& smartPtr);

    /** apply the filter to create a new mesh */
    Mesh apply(const Mesh& inputMesh) const ;

    const bool& serializeLocal() const {return serializeLocal_;}

    bool& serializeLocal() {return serializeLocal_;}
  private:
    bool serializeLocal_;
    
  };
}

#endif
