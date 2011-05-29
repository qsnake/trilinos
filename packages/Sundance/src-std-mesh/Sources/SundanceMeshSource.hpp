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

#ifndef SUNDANCE_MESHSOURCE_H
#define SUNDANCE_MESHSOURCE_H

#include "SundanceDefs.hpp"
#include "SundanceMeshSourceBase.hpp"

#include "SundanceHandle.hpp"

namespace Sundance
{
  /**
   * MeshSource is the user-level interface for objects such as
   * mesh generators and mesh file readers. A MeshSource can
   * create a mesh object with the <tt> getMesh() </tt> method,
   * and if node and element attributes are available, it can
   * access them with the getAttributes() method.
   *
   * Example: read input from a file celled "meshFile" 
   * in Shewchuk's Triangle format, and
   * create a mesh of type BasicSimplicialMesh distributed over the
   * MPI communicator MPI_COMM_WORLD.
   * \code
   * MeshType meshType = new BasicSimplicialMeshType();
   * MeshSource meshSrc = new TriangleMeshReader("meshFile", meshType, MPIComm::world());
   * \endcode
   */
  class MeshSource : public Sundance::Handle<MeshSourceBase>
  {
  public:
    /** Construct an empty mesh source object */
    MeshSource();

    /** Construct from a raw pointer to a mesh source subtype */
    MeshSource(Sundance::Handleable<MeshSourceBase>* rawPtr);

    /** Construct from a smart pointer to a mesh source subtype */
    MeshSource(const RCP<MeshSourceBase>& smartPtr);

    /** Create and return a mesh */
    Mesh getMesh() const ;

    /** Get any attributes associated with the nodes and elements in the
     * mesh. If no attributes exist, the arrays are empty. If the mesh
     * does not exist, it will be created with a cell to getMesh(). */
    void getAttributes(RCP<Array<Array<double> > >& nodeAttributes,
                       RCP<Array<Array<double> > >& elemAttributes) const ;

    /** Return the mesh type to be used by default if no MeshType
     * is given in a MeshSource subtype ctor. The default mesh type
     * can be set by including a specifer such as
     * <pre>
     * <DefaultMesh type="BasicSimplicial"/>
     * </pre>
     * as a child in the XML configuration file. 
     */
    static MeshType& defaultMeshType() ;

    /** access to the MPI communicator */
    const MPIComm& comm() const ;

    static bool& staggerOutput() {static bool rtn=false; return rtn;}

  private:
  };
}

#endif
