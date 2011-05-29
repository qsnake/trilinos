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

#ifndef SUNDANCE_PARTITIONEDLINEMESHER_H
#define SUNDANCE_PARTITIONEDLINEMESHER_H

#include "SundanceDefs.hpp"
#include "SundanceMeshSourceBase.hpp"

namespace Sundance
{
using namespace Teuchos;
  
/**
 * PartitionedLineMesher meshes the one-dimensional interval 
 * \f$\left[a_x, b_x\right]\f$
 * with \f$n_x\f$ elements per processor. 
 */
class PartitionedLineMesher : public MeshSourceBase
{
public:
  /** 
   * Set up a mesher for the interval \f$\left[a_x, b_x\right]\f$
   * with \f$n_x\f$ elements per processor. 
   */
  PartitionedLineMesher(double ax, double bx, int nx,
    const MeshType& meshType,
    const MPIComm& comm = MPIComm::world())
    : 
    MeshSourceBase(meshType, comm),
    ax_(ax), bx_(bx), nx_(nx) {;}

  /** Create a line mesher from a ParameterList */
  PartitionedLineMesher(const ParameterList& params);
    
  /** */
  virtual ~PartitionedLineMesher() {;}

  /** Print a short descriptive std::string */
  virtual std::string description() const 
    {return "PartitionedLineMesher[ax=" + Teuchos::toString(ax_)
        + ", bx=" + Teuchos::toString(bx_)
        + ", nx=" + Teuchos::toString(nx_) + "]";}
      

  /** Return a ref count pointer to self */
  virtual RCP<MeshSourceBase> getRcp() {return rcp(this);}

private:

  /** */
  virtual Mesh fillMesh() const ;

  /** */
  double ax_;
  /** */
  double bx_;
  /** */
  int nx_;

};
}

#endif
