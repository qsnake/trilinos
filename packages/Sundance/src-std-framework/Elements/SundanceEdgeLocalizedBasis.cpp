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

#include "SundanceEdgeLocalizedBasis.hpp"
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


EdgeLocalizedBasis::EdgeLocalizedBasis()
{}

bool EdgeLocalizedBasis::supportsCellTypePair(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  switch(maximalCellType)
  {
    case TriangleCell:
    case TetCell:
    case QuadCell:
    case BrickCell:
      switch(cellType)
      {
        case LineCell:
          return true;
        default:
          return false;
      }
    default:
      return false;
  }
}

void EdgeLocalizedBasis::print(std::ostream& os) const 
{
  os << "EdgeLocalizedBasis()";
}

int EdgeLocalizedBasis::nReferenceDOFsWithoutFacets(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  switch(cellType)
  {
    case LineCell:
      return 1;
    default:
      return 0;
  }
}

void EdgeLocalizedBasis::getReferenceDOFs(
  const CellType& maximalCellType,
  const CellType& cellType,
  Array<Array<Array<int> > >& dofs) const 
{
  typedef Array<int> Aint;
  switch(cellType)
  {
    case LineCell:
      dofs.resize(2);
      dofs[0] = Array<Array<int> >();
      dofs[1] = tuple<Aint>(tuple<int>(0));
      return;
    case TriangleCell:
      dofs.resize(3);
      dofs[0] = tuple(Array<int>());
      dofs[1] = tuple<Aint>(tuple(0), tuple(1), tuple(2));
      dofs[2] = tuple(Array<int>());
      return;
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError, "Cell type "
        << cellType << " not implemented in EdgeLocalizedBasis basis");
  }
}



void EdgeLocalizedBasis::refEval(
  const CellType& cellType,
  const Array<Point>& pts,
  const SpatialDerivSpecifier& sds,
  Array<Array<Array<double> > >& result,
  int verbosity) const
{
  TEST_FOR_EXCEPTION(!(sds.isPartial() || sds.isIdentity()), 
    RuntimeError,
    "cannot evaluate spatial derivative " << sds << " on EdgeLocalizedBasis basis");
  const MultiIndex& deriv = sds.mi();
  typedef Array<double> Adouble;
  result.resize(1);
  result[0].resize(pts.length());

  int dim = dimension(cellType);
  
  if (dim==0)
  {
    result[0] = tuple<Adouble>(tuple(1.0));
  }
  else if (dim==1)
  {
    for (int i=0; i<pts.length(); i++)
    {
      evalOnLine(pts[i], deriv, result[0][i]);
    }
  }
  else if (dim==2)
  {
    for (int i=0; i<pts.length(); i++)
    {
      evalOnTriangle(pts[i], deriv, result[0][i]);
    }
  }
  else if (dim==3)
  {
    for (int i=0; i<pts.length(); i++)
    {
      evalOnTet(pts[i], deriv, result[0][i]);
    }
  }
}

/* ---------- evaluation on different cell types -------------- */


void EdgeLocalizedBasis::evalOnLine(const Point& pt, 
													const MultiIndex& deriv,
													Array<double>& result) const
{
	ADReal one(1.0, 1);
	result.resize(1);
	Array<ADReal> tmp(result.length());

  tmp[0] = one;

	for (int i=0; i<tmp.length(); i++)
		{
			if (deriv.order()==0) result[i] = tmp[i].value();
			else result[i] = tmp[i].gradient()[deriv.firstOrderDirection()];
		}
}

void EdgeLocalizedBasis::evalOnTriangle(const Point& pt, 
  const MultiIndex& deriv,
  Array<double>& result) const
{
  ADReal x = ADReal(pt[0], 0, 2);
	ADReal y = ADReal(pt[1], 1, 2);
	ADReal one(1.0, 2);
	ADReal zero(0.0, 2);

  Array<ADReal> tmp;

  SUNDANCE_OUT(this->verb() > 3, "x=" << x.value() << " y="
    << y.value());

  result.resize(3);
  tmp.resize(3);

  bool onEdge0 = std::fabs(pt[1]) < 1.0e-14;
  bool onEdge1 = std::fabs(1.0-pt[0]-pt[1]) < 1.0e-14;
  bool onEdge2 = std::fabs(pt[0]) < 1.0e-14;
  
  TEST_FOR_EXCEPTION(!(onEdge0 || onEdge1 || onEdge2),
    RuntimeError,
    "EdgeLocalizedBasis should not be evaluated at points not on edges");
  
  TEST_FOR_EXCEPTION((onEdge0 && onEdge1) || (onEdge1 && onEdge2)
    || (onEdge2 && onEdge0), RuntimeError,
    "Ambiguous edge in EdgeLocalizedBasis::evalOnTriangle()");

  if (onEdge0)
  {
    tmp[0] = one;
    tmp[1] = zero;
    tmp[2] = zero;
  }
  if (onEdge1)
  {
    tmp[0] = zero;
    tmp[1] = one;
    tmp[2] = zero;
  }
  if (onEdge2)
  {
    tmp[0] = zero;
    tmp[1] = zero;
    tmp[2] = one;
  }


	for (int i=0; i<tmp.length(); i++)
  {
    SUNDANCE_OUT(this->verb() > 3,
      "tmp[" << i << "]=" << tmp[i].value() 
      << " grad=" << tmp[i].gradient());
    if (deriv.order()==0) result[i] = tmp[i].value();
    else 
      result[i] = tmp[i].gradient()[deriv.firstOrderDirection()];
  }
  
}




void EdgeLocalizedBasis::evalOnTet(const Point& pt, 
  const MultiIndex& deriv,
  Array<double>& result) const
{
  TEST_FOR_EXCEPTION(true, RuntimeError,
    "EdgeLocalizedBasis::evalOnTet not implemented");
}
