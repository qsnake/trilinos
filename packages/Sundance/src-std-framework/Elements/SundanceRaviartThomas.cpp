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

#include "SundanceRaviartThomas.hpp"
#include "SundancePoint.hpp"
#include "SundanceCellType.hpp"
#include "SundanceSpatialDerivSpecifier.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceTypeUtils.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceOut.hpp"
#include "SundanceADReal.hpp"

using namespace Sundance;
using namespace Teuchos;


RaviartThomas::RaviartThomas(int spatialDim)
  : HDivVectorBasis(spatialDim)
{}

bool RaviartThomas::supportsCellTypePair(
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

std::string RaviartThomas::description() const 
{
  return "RaviartThomas()";
}

int RaviartThomas::nReferenceDOFsWithoutFacets(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  if (dimension(cellType) == dimension(maximalCellType)-1)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

void RaviartThomas::getReferenceDOFs(
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
        << cellType << " not implemented in RaviartThomas basis");
  }
}



bool RaviartThomas::lessThan(const BasisDOFTopologyBase* other) const 
{
  if (typeLessThan(this, other)) return true;
  if (typeLessThan(other, this)) return false;

  return false;
}


void RaviartThomas::refEval(
  const CellType& cellType,
  const Array<Point>& pts,
  const SpatialDerivSpecifier& sds,
  Array<Array<Array<double> > >& result,
  int verbosity) const
{
  const MultiIndex& deriv = sds.mi();
  switch(cellType) {
    case PointCell:
    {
      TEST_FOR_EXCEPTION(true, RuntimeError, "evaluation of RaviartThomas elements for PointCell not supported");
    }
    break;
    case LineCell:
    {
      result.resize(1);
      result[0].resize(pts.length());
      Array<ADReal> tmp;
      tmp.resize(2);
      for (int i=0;i<pts.length();i++) {
        ADReal x = ADReal( pts[i][0] , 0 , 1 );
        ADReal one(1.0,1);
        result[0][i].resize(2);
        tmp[0] = one - x;
        tmp[1] = x;
        if (deriv.order()==0) {
          for (int j=0;j<tmp.length();j++) {
            result[0][i][j] = tmp[j].value();
          }
        }
        else {
          for (int j=0;j<tmp.length();j++) {
            result[0][i][j] = tmp[j].gradient()[deriv.firstOrderDirection()];
          }
        }
      }
    }
    break;
    case TriangleCell:
    {
      result.resize(2);
      result[0].resize(pts.length());
      result[1].resize(pts.length());
      Array<ADReal> tmp0;
      Array<ADReal> tmp1;
      tmp0.resize(3);
      tmp1.resize(3);
      for (int i=0;i<pts.length();i++) {
        ADReal x = ADReal(pts[i][0],0,2);
        ADReal y = ADReal(pts[i][1],1,2);
        ADReal one(1.0,2);
        ADReal rt2(std::sqrt(2.0),2);
        result[0][i].resize(3);
        result[1][i].resize(3);
        tmp0[0] = x;
        tmp1[0] = y - one;
        tmp0[1] = rt2 * x;
        tmp1[1] = rt2 * y;
        tmp0[2] = x - one;
        tmp1[2] = y;
        if (deriv.order()==0) {
          for (int j=0;j<tmp0.length();j++) {
            result[0][i][j] = tmp0[j].value();
            result[1][i][j] = tmp1[j].value();
          }
        }
        else {
          for (int j=0;j<tmp0.length();j++) {
            result[0][i][j] = tmp0[j].gradient()[deriv.firstOrderDirection()];
            result[1][i][j] = tmp1[j].gradient()[deriv.firstOrderDirection()];
          }
        }
      }
    }
    break;
    case TetCell:
    {
      TEST_FOR_EXCEPTION(true, RuntimeError, "evaluation of RaviartThomas elements for TetCell not supported");
    }
    break;
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError, "evaluation of RaviartThomas elements unknown cell type");
      break;
  }

  return;
}


void RaviartThomas::print(std::ostream& os) const 
{
  os << "RaviartThomas(" << dim() << ")";
}
