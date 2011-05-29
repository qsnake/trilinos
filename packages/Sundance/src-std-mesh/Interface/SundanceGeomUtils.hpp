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

#ifndef SUNDANCE_GEOMUTILS_H
#define SUNDANCE_GEOMUTILS_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include <list>

namespace Sundance
{
/** \relates Mesh 
 * Return the LID of the maximal cell that contains the point x.
 * If the point lays on an edge between two or more
 * cells, this will return the first of the those cells encountered
 * on a breadth-first search starting at the initial guess.
 */
int findEnclosingCell(const Mesh& mesh, 
  int cellDim,
  int initialGuessLID, 
  const double* x);

/** \relates Mesh
 * Test whether a point is enclosed in a cell. The cell is assumed
 * to be simplicial.
 */
bool cellContainsPoint(const Mesh& mesh, 
  int cellDim,
  int cellLID, const double* x,
  Array<int>& facetLID);

/** \relates Mesh
 * Tests orientation of a point relative to a line.
 */
double orient2D(const double* a, const double* b, const double* x);

/** \relates Mesh 
 * Get the list of maximal neighbors of a cell 
 */
void maximalNeighbors(const Mesh& mesh, int cellDim,
  int cellLID, const Array<int>& facetLID,
  std::list<int>& rtn);

/** \relates Mesh 
 * Pullback a point to local coordinates within a cell
 */
Point pullback(const Mesh& mesh, int cellDim, int cellLID, const double* x);

/** */
void printCell(const Mesh& mesh, int cellLID);

/** */
double volume(const Mesh& mesh, int cellDim, int cellLID);
  

}


#endif


