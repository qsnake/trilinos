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

#include "SundanceBasisDOFTopologyBase.hpp"
#include "SundanceExceptions.hpp"


using namespace Sundance;
using namespace Teuchos;


int BasisDOFTopologyBase::nReferenceDOFsWithFacets(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const 
{
  switch(cellType)
  {
    case TetCell:
      return nReferenceDOFsWithoutFacets(maximalCellType, TetCell)
        +  4*nReferenceDOFsWithoutFacets(maximalCellType, TriangleCell)
        +  6*nReferenceDOFsWithoutFacets(maximalCellType, LineCell)
        +  4*nReferenceDOFsWithoutFacets(maximalCellType, PointCell);
    case BrickCell:
      return nReferenceDOFsWithoutFacets(maximalCellType, BrickCell)
        +  6*nReferenceDOFsWithoutFacets(maximalCellType, QuadCell)
        +  12*nReferenceDOFsWithoutFacets(maximalCellType, LineCell)
        +  8*nReferenceDOFsWithoutFacets(maximalCellType, PointCell);
    case TriangleCell:
      return nReferenceDOFsWithoutFacets(maximalCellType, TriangleCell)
        +  3*nReferenceDOFsWithoutFacets(maximalCellType, LineCell)
        +  3*nReferenceDOFsWithoutFacets(maximalCellType, PointCell);
    case QuadCell:
      return nReferenceDOFsWithoutFacets(maximalCellType, QuadCell)
        +  4*nReferenceDOFsWithoutFacets(maximalCellType, LineCell)
        +  4*nReferenceDOFsWithoutFacets(maximalCellType, PointCell);
    case LineCell:
      return nReferenceDOFsWithoutFacets(maximalCellType, LineCell)
        +  2*nReferenceDOFsWithoutFacets(maximalCellType, PointCell);
    case PointCell:
      return nReferenceDOFsWithoutFacets(maximalCellType, PointCell);
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError,
        "case cellType=" << cellType << " not defined");
  }
  return -1; // -Wall
}
