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

#ifndef SUNDANCE_PERIODIC_LINE_MESHER_HPP
#define SUNDANCE_PERIODIC_LINE_MESHER_HPP

#include "SundanceDefs.hpp"
#include "SundanceMeshReaderBase.hpp"
#include "SundancePeriodicMesh1D.hpp"

namespace Sundance
{
using namespace Teuchos;
  
/**
 * Create a mesh having one triangle
 */
class PeriodicLineMesher : public MeshSourceBase
{
public:
  /** */
  PeriodicLineMesher(
    const double& a, 
    const double& b,
    int nx,
    const MeshType& meshType)
    : MeshSourceBase(meshType, MPIComm::self()),
      nx_(nx),
      a_(a),
      b_(b)
    {}

  /** virtual dtor */
  virtual ~PeriodicLineMesher(){;}


  /** Create a mesh */
  virtual Mesh fillMesh() const 
    {
      RCP<MeshBase> rtn = rcp(new PeriodicMesh1D(a_, b_, nx_));
      return rtn;
    }

  /** Print a short descriptive std::string */
  virtual std::string description() const 
    {return "PeriodicLineMesher";}
      

  /** Return a ref count pointer to self */
  virtual RCP<MeshSourceBase> getRcp() {return rcp(this);}

private:
  int nx_;
  double a_;
  double b_;
};

}

#endif
