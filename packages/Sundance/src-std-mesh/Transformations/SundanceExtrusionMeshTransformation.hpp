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

#ifndef SUNDANCE_EXTRUSIONMESHFILTER_H
#define SUNDANCE_EXTRUSIONMESHFILTER_H

#include "SundanceDefs.hpp"
#include "SundanceMeshTransformationBase.hpp"

namespace Sundance
{
  
  /**
   * ExtrusionMeshTransformation extrudes a 2D mesh to 3D. 
   */
  class ExtrusionMeshTransformation : public MeshTransformationBase
  {
  public:
    /** Construct a filter to extrude a 2D mesh from
     * the plane \f$z=z_0\f$ to the plane \f$z=z_1\f$ in 
     * \f$n_z\f$ steps. */
    ExtrusionMeshTransformation(double z0, double z1, int nzLevels,
                        const MeshType& meshType)
      : MeshTransformationBase(meshType),
        z0_(z0), z1_(z1), nzLevels_(nzLevels) {;}

    /** virtual dtor */
    virtual ~ExtrusionMeshTransformation() {;}

    /** Apply the filter to an input mesh, returning an output mesh */
    virtual Mesh apply(const Mesh& inputMesh) const ;

    /** Print a short descriptive std::string */
    virtual std::string description() const 
    {return "ExtrusionMeshTransformation[z0=" + Teuchos::toString(z0_)
       + ", z1=" + Teuchos::toString(z1_)
       + ", nz=" + Teuchos::toString(nzLevels_) + "]";}
      

    /** Return a ref count pointer to self */
    virtual RCP<MeshTransformationBase> getRcp() {return rcp(this);}

  private:
    
    /** */
    double z0_;
    
    /** */
    double z1_;

    /** */
    int nzLevels_;
    
  };
}

#endif
