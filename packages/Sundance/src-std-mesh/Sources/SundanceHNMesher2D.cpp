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


#include "SundanceHNMesher2D.hpp"

using namespace Sundance;
using namespace Teuchos;

REFINE_MESH_ESTIMATE( DummyNoRefineClass , { return false; } , {return 0;} )
MESH_DOMAIN( DummyAllMeshDomain , {return true;})

const RefinementClass HNMesher2D::dummyRefineClass_ = (new DummyNoRefineClass());

const MeshDomainDef HNMesher2D::dummyMeshDomain_ = (new DummyAllMeshDomain());

HNMesher2D::HNMesher2D(const ParameterList& params)
: MeshSourceBase(params),
  _position_x(params.get<double>("position_x")),
  _position_y(params.get<double>("position_y")),
  _offset_x(params.get<double>("offset_x")),
  _offset_y(params.get<double>("offset_y")),
  _resolution_x(params.get<int>("resolution_x")),
  _resolution_y(params.get<int>("resolution_y")),
  refineClass_(HNMesher2D::dummyRefineClass_),
  meshDomain_(HNMesher2D::dummyMeshDomain_)
{
	// nothing to do here
}

Mesh HNMesher2D::fillMesh() const
{
	// here we create the mesh and return to the Sundance
	HNMesh2D *hnodegrid;
	Mesh mesh = createMesh(2);
	// get the pointer
	hnodegrid = (HNMesh2D*) mesh.ptr().get();
	hnodegrid->createMesh( _position_x , _position_y , _offset_x , _offset_y , _resolution_x , _resolution_y ,
			                refineClass_ , meshDomain_);
	return (mesh);
}

