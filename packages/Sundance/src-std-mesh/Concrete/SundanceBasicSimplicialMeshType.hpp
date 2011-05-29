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

#ifndef SUNDANCE_BASICSIMPLICIALMESHTYPE_H
#define SUNDANCE_BASICSIMPLICIALMESHTYPE_H


#include "SundanceDefs.hpp"
#include "SundanceMeshBase.hpp"
#include "SundanceBasicSimplicialMesh.hpp"

namespace Sundance
{
using namespace Teuchos;
  
/**
 * BasicSimplicialMeshType is used to create
 * BasicSimplicialMesh objects.
 */
class BasicSimplicialMeshType : public MeshTypeBase
{
public:
  /** Empty ctor */
  BasicSimplicialMeshType(const MeshEntityOrder& order=ExodusMeshOrder) 
    : order_(order) {;}

  /** virtual dtor */
  virtual ~BasicSimplicialMeshType(){;}

  /** Create a mesh of the given dimension */
  virtual RCP<MeshBase> createEmptyMesh(int dim,
    const MPIComm& comm) const 
    {return rcp(new BasicSimplicialMesh(dim, comm, order_));}

  /** */
  std::string description() const {return "BasicSimplicialMeshType";}

  /** Return a ref count pointer to self */
  virtual RCP<MeshTypeBase> getRcp() {return rcp(this);}

private:
  MeshEntityOrder order_;
      
};
}
#endif
