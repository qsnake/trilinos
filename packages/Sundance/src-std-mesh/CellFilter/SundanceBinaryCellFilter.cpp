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

#include "SundanceBinaryCellFilter.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOrderedTuple.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

BinaryCellFilter::BinaryCellFilter(const CellFilter& left,
                                   const CellFilter& right,
                                   const CellFilterOpType& op)
  : CellFilterBase(), op_(op), left_(left), right_(right)
{
  std::string str;
  switch(op)
  {
    case Union:
      str = "Union(";
      break;
    case Intersection:
      str = "Intersection(";
      break;
    default:
      str = "SetDifference(";
    }

  setName(str + left.toString() + ", " + right.toString() + ")");
}

int BinaryCellFilter::dimension(const Mesh& mesh) const
{
  int d1 = left_.dimension(mesh);
  int d2 = right_.dimension(mesh);

  TEST_FOR_EXCEPTION(d1 != d2, RuntimeError,
                     "BinaryCellFilter::dimension() mismatched dimensions. "
                     "Left filter has dimension d1=" << d1 << " but "
                     "right filter has dimension d2=" << d2);

  return d1;
}

CellSet BinaryCellFilter::internalGetCells(const Mesh& mesh) const
{
  Tabs tab;

  SUNDANCE_OUT(this->verb() > 2,
               "cell filter " << toXML().toString() << " is getting its cells for mesh " 
               << mesh.id());

  CellSet L = left_.getCells(mesh);
  CellSet R = right_.getCells(mesh);


  SUNDANCE_OUT(this->verb() > 2,
               "cell filter " << toXML().toString() << " is performing its operation");
  
  switch(op_)
    {
    case Union:
      return L.setUnion(R);
    case Intersection:
      return L.setIntersection(R);
    case Difference:
      return L.setDifference(R);
    }
  TEST_FOR_EXCEPTION(true, RuntimeError, "unknown cell filter op type" << op_ 
                     << " in BinaryCellFilter::internalGetCells()");
  return L; // -Wall
}

string BinaryCellFilter::opName() const 
{
  switch(op_)
    {
    case Union:
      return "UnionCellFilter";
    case Intersection:
      return "IntersectionCellFilter";
    case Difference:
      return "DifferenceCellFilter";
    }
  TEST_FOR_EXCEPTION(true, RuntimeError, "unknown cell filter op type" << op_ 
                     << " in BinaryCellFilter::opName()");
  return "UnionCellFilter";//-Wall
}

XMLObject BinaryCellFilter::toXML() const 
{
  XMLObject rtn(opName());
  rtn.addChild(left_.toXML());
  rtn.addChild(right_.toXML());
  return rtn;
}

bool BinaryCellFilter::lessThan(const CellFilterStub* other) const
{
  const BinaryCellFilter* B 
    = dynamic_cast<const BinaryCellFilter*>(other);

  TEST_FOR_EXCEPTION(B==0,
                     InternalError,
                     "argument " << other->toXML() 
                     << " to BinaryCellFilter::lessThan() should be "
                     "a BinaryCellFilter pointer.");

  return OrderedTriple<CellFilterOpType, CellFilter, CellFilter>(op_, left_, right_) 
    < OrderedTriple<CellFilterOpType, CellFilter, CellFilter>(B->op_, B->left_, B->right_) ;
}
