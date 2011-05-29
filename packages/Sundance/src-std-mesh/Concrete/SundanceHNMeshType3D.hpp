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

/*
 * SundanceHNMeshType3D.hpp
 *
 *  Created on: May 30, 2010
 *      Author: benk
 */

#ifndef SUNDANCE_HNMESHTYPE3D_HPP_
#define SUNDANCE_HNMESHTYPE3D_HPP_

#include "SundanceDefs.hpp"
#include "SundanceMeshBase.hpp"
#include "SundanceHNMesh3D.hpp"

namespace Sundance
{
  using namespace Teuchos;

  /**
   * Class for a simple capable to refine, and to work in parallel structured 3D mesh type <br>
   * Mesh is based on the trisection refinement.
   */
  class HNMeshType3D : public MeshTypeBase
  {
  public:
    /** Empty ctor */
	  HNMeshType3D(const MeshEntityOrder& order=ExodusMeshOrder)
	    : order_(order) {;}

    /** virtual dtor */
    virtual ~HNMeshType3D(){;}

    /** Create a mesh of the given dimension */
    virtual RCP<MeshBase> createEmptyMesh(int dim,
                                          const MPIComm& comm) const
    // this line is never used since we create directly the mesh at the Mesher
    {return rcp(new HNMesh3D(dim, comm, order_));}

    /** */
    std::string description() const {return "HNMeshType3D";}

    /** Return a ref count pointer to self */
    virtual RCP<MeshTypeBase> getRcp() {return rcp(this);}

  private:
    MeshEntityOrder order_;

  };
}

#endif /* SUNDANCEHNMESHTYPE3D_HPP_ */
