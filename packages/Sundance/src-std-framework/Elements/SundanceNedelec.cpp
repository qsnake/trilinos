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

#include "SundanceNedelec.hpp"
#include "SundanceADReal.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceSpatialDerivSpecifier.hpp"
#include "SundancePoint.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;


Nedelec::Nedelec(int spatialDim)
  : HCurlVectorBasis(spatialDim)
{}

bool Nedelec::supportsCellTypePair(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  switch(maximalCellType)
  {
    case TriangleCell:
      switch(cellType)
      {
        case TriangleCell:
        case LineCell:
        case PointCell:
          return true;
        default:
          return false;
      }
    case TetCell:
      // tets not yet implemented
      return false;
      switch(cellType)
      {
        case TetCell:
        case TriangleCell:
        case LineCell:
        case PointCell:
          return true;
        default:
          return false;
      }
    default:
      return false;
  }
}

void Nedelec::print(std::ostream& os) const 
{
  os << "Nedelec()";
}

int Nedelec::nReferenceDOFsWithoutFacets(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  TEST_FOR_EXCEPT(maximalCellType != TriangleCell);
  switch(cellType)
    {
    case PointCell:
      return 0;
    case LineCell:
      return 1;
    case TriangleCell:
      return 0;
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError, "Cell type "
                         << cellType << " not implemented in Nedelec basis");
      return -1; // -Wall
    }
}

void Nedelec::getReferenceDOFs(
  const CellType& maximalCellType,
  const CellType& cellType,
  Array<Array<Array<int> > >& dofs) const 
{
  switch(cellType)
    {
    case PointCell:
      dofs.resize(1);
      dofs[0] = tuple(Array<int>());
      return;
    case LineCell:
      dofs.resize(2);
      dofs[0] = tuple(Array<int>());
      dofs[1] = tuple<Array<int> >(tuple(0));
      return;
    case TriangleCell:
      dofs.resize(3);
      dofs[0] = tuple(Array<int>());
      dofs[1] = tuple<Array<int> >(tuple(0), tuple(1), tuple(2));
      dofs[2] = tuple(Array<int>());
      return;
    case TetCell:
      dofs.resize(4);
      dofs[0] = tuple(Array<int>());
      dofs[1] = tuple<Array<int> >(tuple(0), tuple(1), tuple(2), tuple(3));
      dofs[2] = tuple(Array<int>());
      dofs[3] = tuple(Array<int>());
      return;
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError, "Cell type "
                         << cellType << " not implemented in Nedelec basis");
    }
}



void Nedelec::refEval(
  const CellType& cellType,
  const Array<Point>& pts,
  const SpatialDerivSpecifier& deriv,
  Array<Array<Array<double> > >& result,
  int verbosity) const
{
  TEST_FOR_EXCEPTION(true, RuntimeError, "evaluation of Nedelec elements not yet supported");
}


