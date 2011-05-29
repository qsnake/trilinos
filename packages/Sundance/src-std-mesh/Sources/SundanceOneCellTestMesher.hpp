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

#ifndef SUNDANCE_ONE_CELL_TEST_MESHER_H
#define SUNDANCE_ONE_CELL_TEST_MESHER_H

#include "SundanceDefs.hpp"
#include "SundanceMeshReaderBase.hpp"

namespace Sundance
{
using namespace Teuchos;
  
/**
 * Create a mesh having one triangle
 */
class OneTriangleMesher : public MeshSourceBase
{
public:
  /** */
  OneTriangleMesher(
    const Point& A, 
    const Point& B,
    const Point& C,
    const MeshType& meshType);

  /** virtual dtor */
  virtual ~OneTriangleMesher(){;}


  /** Create a mesh */
  virtual Mesh fillMesh() const ;

  /** Print a short descriptive std::string */
  virtual std::string description() const 
    {return "OneTriangleMesher";}
      

  /** Return a ref count pointer to self */
  virtual RCP<MeshSourceBase> getRcp() {return rcp(this);}

private:
  Point A_;
  Point B_;
  Point C_;
};

/**
 * Create a mesh having one tet
 */
class OneTetMesher : public MeshSourceBase
{
public:
  /** */
  OneTetMesher(
    const Point& A, 
    const Point& B,
    const Point& C,
    const Point& D,
    const MeshType& meshType);

  /** virtual dtor */
  virtual ~OneTetMesher(){;}


  /** Create a mesh */
  virtual Mesh fillMesh() const ;

  /** Print a short descriptive std::string */
  virtual std::string description() const 
    {return "OneTetMesher";}
      

  /** Return a ref count pointer to self */
  virtual RCP<MeshSourceBase> getRcp() {return rcp(this);}

private:
  Point A_;
  Point B_;
  Point C_;
  Point D_;
};
}

#endif
