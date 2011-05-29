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


#include "SundanceOneCellTestMesher.hpp"

using namespace Sundance;
using namespace Teuchos;


OneTriangleMesher::OneTriangleMesher(
  const Point& A,
  const Point& B,
  const Point& C,
  const MeshType& meshType)
  : MeshSourceBase(meshType, MPIComm::world()),
    A_(A), B_(B), C_(C)
{}


Mesh OneTriangleMesher::fillMesh() const
{
  Mesh mesh = createMesh(2);
  
  int a = mesh.addVertex(0, A_, 0, 0);
  int b = mesh.addVertex(1, B_, 0, 0);
  int c = mesh.addVertex(2, C_, 0, 0);

  mesh.addElement(0, tuple(a,b,c), 0, 0);

  return mesh;
}
  


OneTetMesher::OneTetMesher(
  const Point& A,
  const Point& B,
  const Point& C,
  const Point& D,
  const MeshType& meshType)
  : MeshSourceBase(meshType, MPIComm::world()),
    A_(A), B_(B), C_(C), D_(D)
{}


Mesh OneTetMesher::fillMesh() const
{
  Mesh mesh = createMesh(3);
  
  int a = mesh.addVertex(0, A_, 0, 0);
  int b = mesh.addVertex(1, B_, 0, 0);
  int c = mesh.addVertex(2, C_, 0, 0);
  int d = mesh.addVertex(3, D_, 0, 0);

  mesh.addElement(0, tuple(a,b,c,d), 0, 0);

  return mesh;
}
  
