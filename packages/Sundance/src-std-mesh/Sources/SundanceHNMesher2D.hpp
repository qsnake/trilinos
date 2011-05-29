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

#ifndef SUNDANCE_HNMESHER2D_H_
#define SUNDANCE_HNMESHER2D_H_

#include "SundanceDefs.hpp"
#include "SundanceMeshSourceBase.hpp"
#include "SundanceHNMesh2D.hpp"
#include "SundanceRefinementBase.hpp"
#include "SundanceRefinementClass.hpp"
#include "SundanceDomainDefinition.hpp"

namespace Sundance
{
  /** forward declaration */
  class RefinementClass;
  class MeshDomainDef;

class HNMesher2D : public MeshSourceBase
  {
  public:
    /**     */
	HNMesher2D(
			   double position_x, double position_y,
               double offset_x , double offset_y,
               int resolution_x , int resolution_y ,
               const MeshType& meshType,
               const MPIComm& comm = MPIComm::world())
      :
      MeshSourceBase(meshType, comm),
      _position_x(position_x), _position_y(position_y),
      _offset_x(offset_x), _offset_y(offset_y),
      _resolution_x(resolution_x) , _resolution_y(resolution_y) ,
    	refineClass_(dummyRefineClass_) ,
    	meshDomain_(dummyMeshDomain_) {;}

    /**     */
	HNMesher2D(
			   double position_x, double position_y,
               double offset_x , double offset_y,
               int resolution_x , int resolution_y ,
               const MeshType& meshType,
               const RefinementClass& refineClass ,
               const MPIComm& comm = MPIComm::world())
      :
      MeshSourceBase(meshType, comm),
      _position_x(position_x), _position_y(position_y),
      _offset_x(offset_x), _offset_y(offset_y),
      _resolution_x(resolution_x) , _resolution_y(resolution_y) ,
    	refineClass_(refineClass) ,
    	meshDomain_(dummyMeshDomain_) {;}

    /**     */
	HNMesher2D(
			   double position_x, double position_y,
               double offset_x , double offset_y,
               int resolution_x , int resolution_y ,
               const MeshType& meshType,
               const MeshDomainDef& meshDomain ,
               const MPIComm& comm = MPIComm::world())
      :
      MeshSourceBase(meshType, comm),
      _position_x(position_x), _position_y(position_y),
      _offset_x(offset_x), _offset_y(offset_y),
      _resolution_x(resolution_x) , _resolution_y(resolution_y) ,
    	refineClass_(dummyRefineClass_) ,
    	meshDomain_(meshDomain) {;}

    /**     */
	HNMesher2D(
			   double position_x, double position_y,
               double offset_x , double offset_y,
               int resolution_x , int resolution_y ,
               const MeshType& meshType,
               const RefinementClass& refineClass ,
               const MeshDomainDef& meshDomain ,
               const MPIComm& comm = MPIComm::world())
      :
      MeshSourceBase(meshType, comm),
      _position_x(position_x), _position_y(position_y),
      _offset_x(offset_x), _offset_y(offset_y),
      _resolution_x(resolution_x) , _resolution_y(resolution_y) ,
    	refineClass_(refineClass) ,
    	meshDomain_(meshDomain) {;}

    /** Create a rectangle mesher from a ParameterList */
	HNMesher2D(const ParameterList& params);

    /** */
    virtual ~HNMesher2D() {;}

    /** Print a short descriptive std::string */
    virtual std::string description() const
    {return "HNMesher2D[pos x =" + Teuchos::toString(_position_x)
       + ", pos y=" + Teuchos::toString(_position_y)
       + ", offset x=" + Teuchos::toString(_offset_x) +
       + ", offset y=" + Teuchos::toString(_offset_y)
       + ", resolution_x=" + Teuchos::toString(_resolution_x)
       + ", resolution_y=" + Teuchos::toString(_resolution_y)+"]";}


    /** Return a ref count pointer to self */
    virtual RCP<MeshSourceBase> getRcp() {return rcp(this);}


  protected:

    /** The method which all Mesher should have */
    virtual Mesh fillMesh() const ;

  private:

    /** X coordinate of the origin point (lower left)*/
    double _position_x;
    /** Y coordinate of the origin point (lower left)*/
    double _position_y;
    /** offset (length) of the grid in the X direction*/
    double _offset_x;
    /** offset (length) of the grid in the Y direction*/
    double _offset_y;
    /** On the coarse level the resolution on the X axis */
    int _resolution_x;
    /** On the coarse level the resolution on the Y axis */
    int _resolution_y;

    /** refinement class */
    const RefinementClass refineClass_;

    /** mesh domain (which must not coincide with the whole mesh)*/
    const MeshDomainDef meshDomain_;


    /** static dummy class if the user does not provide refinement class */
    static const RefinementClass dummyRefineClass_;

    /** static domain class if the user does not provide one */
    static const MeshDomainDef dummyMeshDomain_;
  };
}
#endif /* SUNDANCE_HNMESHER2D_H_ */
