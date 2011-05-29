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

#include "SundanceCellType.hpp"
#include "SundanceExceptions.hpp"

using namespace Sundance;
using namespace Sundance;

namespace Sundance
{

  std::string toString(const CellType& cellType)
  {
    switch(cellType)
      {
      case NullCell:
        return "NullCell";
      case PointCell:
        return "PointCell";
      case LineCell:
        return "LineCell";
      case TriangleCell:
        return "TriangleCell";
      case QuadCell:
        return "QuadCell";
      case TetCell:
        return "TetCell";
      case BrickCell:
        return "BrickCell";
      case PrismCell:
        return "PrismCell";
      }
    return "NullCell"; // -Wall
  }

  int dimension(const CellType& cellType)
  {
    switch(cellType)
      {
      case NullCell:
        return -1;
      case PointCell:
        return 0;
      case LineCell:
        return 1;
      case TriangleCell:
      case QuadCell:
        return 2;
      case TetCell:
      case BrickCell:
      case PrismCell:
        return 3;
      }
    return -1; // -Wall
  }

  int numFacets(const CellType& cellType, int facetDim)
  {
    int d = dimension(cellType);
    if (facetDim == d) return 1;

    TEST_FOR_EXCEPTION(facetDim > d, RuntimeError,
                       "invalid facet dim " << facetDim << " for cell "
                       << toString(cellType));

    switch(cellType)
      {
      case NullCell:
      case PointCell:
        return -1;
      case LineCell:
        return 2;
      case TriangleCell:
        return 3;
      case TetCell:
        if (facetDim==0 || facetDim==2) return 4;
        return 6;
      case QuadCell:
        return 4;
      case BrickCell:
        if (facetDim==0) return 8;
        else if (facetDim==1) return 12;
        else return 6;
      case PrismCell:
        if (facetDim==0) return 6;
        else if (facetDim==1) return 9;
        else return 5;
      }
    return -1; // -Wall
  }


  CellType facetType(const CellType& cellType, int facetDim, int facetIndex)
  {
    int d = dimension(cellType);
    if (facetDim == d) return cellType;
    TEST_FOR_EXCEPTION(facetDim > d, RuntimeError,
                       "invalid facet dim " << facetDim << " for cell "
                       << toString(cellType));

    if (facetDim==0) return PointCell;
    if (facetDim==1) return LineCell;

    switch(cellType)
      {
      case NullCell:
      case PointCell:
      case LineCell:
      case TriangleCell:
      case QuadCell:
        return NullCell;

      case TetCell:
        return TriangleCell;
      case BrickCell:
        return QuadCell;
      case PrismCell:
        if (facetIndex==0 || facetIndex==4) return TriangleCell;
        return QuadCell;
      }
    return NullCell;
  }
  
}
