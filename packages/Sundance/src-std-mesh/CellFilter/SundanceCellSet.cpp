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

#include "SundanceCellSet.hpp"
#include "SundanceExplicitCellSet.hpp"
#include "SundanceImplicitCellSet.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceExceptions.hpp"
#include <algorithm>
#include <iterator>

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;


CellSet::CellSet(const Mesh& mesh, int cellDim,
                 const CellType& cellType,
                 const Set<int>& cellLIDs)
  : Handle<CellSetBase>(rcp(new ExplicitCellSet(mesh, cellDim, cellType, cellLIDs)))
{}

CellSet CellSet::setUnion(const CellSet& other) const
{
  if (isNull()) return other;
  if (other.isNull()) return *this;

  ExplicitCellSet* rtn = new ExplicitCellSet(mesh(), dimension(), cellType());

  checkCompatibility("union", other);
  
  Set<int>& cells = rtn->cells();

  std::set_union(this->begin(), this->end(), other.begin(), other.end(), 
                 std::insert_iterator<Set<int> >(cells, cells.begin()));
  
  return rtn;
}

CellSet CellSet::setIntersection(const CellSet& other) const
{
  if (isNull()) return *this;
  if (other.isNull()) return other;

  ExplicitCellSet* rtn = new ExplicitCellSet(mesh(), dimension(), cellType());

  checkCompatibility("intersection", other);
  
  Set<int>& cells = rtn->cells();

  std::set_intersection(this->begin(), this->end(), other.begin(), other.end(), 
                        std::insert_iterator<Set<int> >(cells, cells.begin()));
  
  return rtn;
}

CellSet CellSet::setDifference(const CellSet& other) const
{
  if (isNull()) return *this;
  if (other.isNull()) return *this;

  ExplicitCellSet* rtn = new ExplicitCellSet(mesh(), dimension(), cellType());

  checkCompatibility("difference", other);
  
  Set<int>& cells = rtn->cells();

  std::set_difference(this->begin(), this->end(), other.begin(), other.end(), 
                      std::insert_iterator<Set<int> >(cells, cells.begin()));
  
  return rtn;
}


void CellSet::checkCompatibility(const std::string& op, const CellSet& other) const 
{
  TEST_FOR_EXCEPTION(meshID() != other.meshID(), RuntimeError,
                     "CellSet::checkCompatibility(): "
                     "incompatible mesh ID numbers in " << op
                     << ". LHS=" << meshID() << " RHS=" << other.meshID());

  TEST_FOR_EXCEPTION(dimension() != other.dimension(), RuntimeError,
                     "CellSet::checkCompatibility() incompatible dimensions in " << op
                     << "LHS has "
                     "dimension=" << dimension() << " but RHS has dimension="
                     << other.dimension());
  
  TEST_FOR_EXCEPTION(cellType() != other.cellType(), RuntimeError,
                     "CellSet::checkCompatibility() incompatible cell types. "
                     " in " << op << " LHS has "
                     "cellType=" << cellType() << " but RHS has cellType="
                     << other.cellType());

  SUNDANCE_OUT(this->verb() > 2,
               "Set operation: " << op << std::endl
               << "LHS cells: " << *this << std::endl
               << "RHS cells: " << other);
               
}


bool CellSet::areFacetsOf(const CellSet& other) const
{
  Array<int> cofacetLIDs;
  int myDim = dimension();
  int cofacetDim = other.dimension();
  CellType cofacetType = other.cellType();
  if (myDim >= cofacetDim) return false;

  for (CellIterator i=begin(); i!=end(); i++)
    {
      int cellLID = *i;
      mesh().getCofacets(myDim, cellLID, cofacetDim, cofacetLIDs);
      Set<int> cofacetSet;
      for (int c=0; c<cofacetLIDs.size(); c++)
        {
          int cf = cofacetLIDs[c];
          cofacetSet.put(cf);
        }
      CellSet myCofacetSet(mesh(), cofacetDim, cofacetType, cofacetSet); 
      CellSet intersection = myCofacetSet.setIntersection(other);      
      /* if the intersection is empty, then we have found a cell
       * that is not a facet of one of the other cells */
      if (intersection.begin()==intersection.end()) return false;
    }
  return true;
}

bool CellSet::operator<(const CellSet& other) const
{
  Tabs tab;
  bool rtn = ptr()->lessThan(other.ptr().get());
  return rtn;
}

