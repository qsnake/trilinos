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

#ifndef SUNDANCE_MESHFILTERBASE_H
#define SUNDANCE_MESHFILTERBASE_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceHandleable.hpp"
#include "Teuchos_Describable.hpp"
#include "SundancePrintable.hpp"
#include "SundanceNoncopyable.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceIncrementallyCreatableMesh.hpp"

namespace Sundance
{
/**
 * MeshSourceBase provides the internal interface for mesh filters, i.e.,
 * objects that take an input mesh and produce a new mesh. Examples
 * of filter operations are refinement, load balancing,
 * and extrusion from 2D to 3D. 
 * The action of a mesh filter should be independent
 * of the internal mesh representation used. To allow user-level
 * specification of the type of internal mesh representation to be
 * used, a MeshTransformationBase is constructed with a MeshType object
 * which acts as a factory to produce empty output meshes.
 *
 * If the
 * communicator has more than one processor, the mesh created will
 * be distributed.
 *
 * <h4> Writing your own MeshTransformationBase subtype </h4>
 *
 * The only method you will need to override is
 * <ul>
 * <li> <tt>virtual Mesh apply(const Mesh& inputMesh) const = 0 </tt> 
 * </ul>
 * which is where you do the filter action and return an output
 * mesh. This method
 * should usually physically create the mesh with a call to createMesh(),
 * ensuring that the correct mesh representation type is created
 * using the MeshType factory with which the filter was constructed.
 *
 * See the ExtrustionMeshTransformation source code for a very simple
 * example of how to write a mesh filter subtype. 
 *
 * Optionally, you can override the description() method to 
 * provide more informative descriptive output than the std::string
 * <tt>"MeshTransformationBase[unknown subtype]".</tt>
 */
class MeshTransformationBase : public Sundance::Handleable<MeshTransformationBase>,
                               public Sundance::Printable,
                               public Teuchos::Describable,
                               public Noncopyable,
                               public ObjectWithClassVerbosity<MeshTransformationBase>
{
public:
  /** Construct with a mesh type, which specifies
   *  the type of mesh to be built when the filter is applied. */
  MeshTransformationBase(const MeshType& meshType)
    : meshType_(meshType) {;}

  /** virtual dtor */
  virtual ~MeshTransformationBase(){;}

      
  /** Apply the filter to the given input mesh, 
   *  producing an output mesh */
  virtual Mesh apply(const Mesh& inputMesh) const = 0 ;

  /** \name Printable interface */
  //@{
  /** Print to a stream */
  virtual void print(std::ostream& os) const {os << description();}
  //@}

  /** \name Describable interface */
  //@{
  /** Print to a stream */
  virtual std::string description() const 
    {return "MeshTransformationBase[unknown subtype]";}
  //@}
      
protected:

  /** createMesh() allocates the mesh object with a call to 
   * meshType's createMesh() method. */
  Mesh createMesh(int dim, const MPIComm& comm) const ;

private:
  /** */
  MeshType meshType_;


};

}


#endif
