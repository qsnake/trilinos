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

#ifndef SUNDANCE_PARTITIONEDRECTANGLEMESHER_H
#define SUNDANCE_PARTITIONEDRECTANGLEMESHER_H

#include "SundanceDefs.hpp"
#include "SundanceMeshSourceBase.hpp"

namespace Sundance
{
using namespace Teuchos;



/**
 * PartitionedRectangleMesher meshes the rectangle 
 * \f$ \left[a_x, b_x\right] \otimes \left[a_y, b_y\right] \f$
 * with \f$ n_x \otimes n_y \f$ elements per processor. The 
 * rectangle is partitioned among processors, with \f$np_x\f$
 * equal sized 
 * subdomains in the \f$x\f$ direction and \f$np_y\f$ in the \f$y\f$
 * direction.
 */
class PartitionedRectangleMesher : public MeshSourceBase
{
public:
  /** 
   * Set up meshing of the rectangle 
   * \f$ \left[a_x, b_x\right] \otimes \left[a_y, b_y\right] \f$
   * with \f$ n_x \otimes n_y \f$ elements per processor. The 
   * rectangle is partitioned among processors, with \f$np_x\f$
   * equal sized 
   * subdomains in the \f$x\f$ direction and \f$np_y\f$ in the \f$y\f$
   * direction.
   */
  PartitionedRectangleMesher(double ax, double bx, int nx, int npx,
    double ay, double by, int ny, int npy,
    const MeshType& meshType,
    const MPIComm& comm = MPIComm::world())
    : 
    MeshSourceBase(meshType, comm),
    ax_(ax), bx_(bx), nx_(nx), npx_(npx),
    ay_(ay), by_(by), ny_(ny), npy_(npy) {;}

    
  /** Create a rectangle mesher from a ParameterList */
  PartitionedRectangleMesher(const ParameterList& params);
    
  /** */
  virtual ~PartitionedRectangleMesher() {;}

  /** Print a short descriptive std::string */
  virtual std::string description() const 
    {return "PartitionedRectangleMesher[ax=" + Teuchos::toString(ax_)
        + ", bx=" + Teuchos::toString(bx_)
        + ", nx=" + Teuchos::toString(nx_) +
        + ", ay=" + Teuchos::toString(ay_)
        + ", by=" + Teuchos::toString(by_)
        + ", ny=" + Teuchos::toString(ny_) + "]";}



  /** Find a nearly equal balance between X and Y partitions */
  static void balanceXY(int n, int* npx, int* npy);
    

  /** Return a ref count pointer to self */
  virtual RCP<MeshSourceBase> getRcp() {return rcp(this);}



protected:

  /** */
  virtual Mesh fillMesh() const ;

private:

  /** */
  double ax_;
  /** */
  double bx_;
  /** */
  int nx_;
  /** */
  int npx_;

  /** */
  double ay_;
  /** */
  double by_;
  /** */
  int ny_;
  /** */
  int npy_;


};
}

#endif
