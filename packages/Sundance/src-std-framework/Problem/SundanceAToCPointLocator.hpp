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

#ifndef SUNDANCE_ATOCPOINTLOCATOR_H
#define SUNDANCE_ATOCPOINTLOCATOR_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceExpr.hpp"
#include "SundanceCellFilter.hpp"

namespace Sundance
{
  using namespace Sundance;
  using namespace Sundance;
  using namespace Sundance;
  using namespace Sundance;
  using namespace Sundance;
  using namespace Teuchos;

  /**
   * AToCPointLocator finds the cell index for a point within an unstructured
   * mesh.
   *
   * Note: not tested in parallel.
   */
  class AToCPointLocator
  {
  public:
    /** */
    AToCPointLocator(const Mesh& mesh, 
                     const CellFilter& subdomain,
                     const std::vector<int>& nx);



    /** Find the index of a point in an overlaid structured grid. */
    int getGridIndex(const double* x) const ;

    /** Use an overlaid structured grid to estimate the location of the point. */
    int guessCell(const double* x) const 
    {return (*table_)[getGridIndex(x)];}

    /** Find the cell that contains the specified point */
    int findEnclosingCell(int initialGuessLID, const double* x) const ;

    /** Find the cell that contains the specified point, also
     * computing local coordinates within that cell. */
    int findEnclosingCell(int initialGuessLID, const double* x,
                          double* localCoords) const ;

    /** */
    void fillMaximalNeighbors(int cellLID, const int* facetLID) const ;

    /** Test whether a point is within a specified cell */
    bool cellContainsPoint(int cellLID, 
                           const double* x, 
                           const int* facetLID) const ;

    /** Test whether a point is within a specified cell, and if so,
     * compute local coordinates within that cell. */
    bool cellContainsPoint(int cellLID, 
                           const double* x, 
                           const int* facetLID,
                           double* localCoords) const ;

    /** */
    const Mesh& mesh() const {return mesh_;}

    /** */
    const CellFilter& subdomain() const {return subdomain_;}



    /** */
    static Point makePoint(int dim, const double* x) ;
  private:

    /** Find the range of structured grid cells within the bounding box of 
     * a cell. */
    void getGridRange(const Mesh& mesh, int cellDim, int cellLID,
                 Array<int>& lowIndex, Array<int>& highIndex) const;




    int dim_;
    Mesh mesh_;
    int nFacets_;
    std::vector<int> nx_;
    Array<double> low_;
    Array<double> high_;
    Array<double> dx_;
    RCP<Array<int> > table_;
    CellFilter subdomain_;
    mutable Array<RCP<Set<int> > > neighborSet_;
  };
}


#endif
