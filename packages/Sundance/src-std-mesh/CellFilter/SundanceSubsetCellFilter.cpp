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

#include "SundanceSubsetCellFilter.hpp"
#include "SundanceExplicitCellSet.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOrderedTuple.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

SubsetCellFilter::SubsetCellFilter(const CellFilter& superset,
                                   const CellPredicate& predicate)
  : CellFilterBase(), superset_(superset), predicate_(predicate)
{
  int verb=0;
  SUNDANCE_MSG3(verb, "creating subset cell filter: [" 
    << predicate.description()
    << "]");
  setName("Subset[sup=" + superset.description() + ", pred="
    + predicate.description()+"]");
}

SubsetCellFilter::~SubsetCellFilter()
{
//  Out::os() << "~SubsetCellFilter()" << std::endl;
}

XMLObject SubsetCellFilter::toXML() const 
{
  XMLObject rtn("SubsetCellFilter");
  rtn.addAttribute("id", Teuchos::toString(id()));
  rtn.addChild(predicate_.toXML());
  return rtn;
}

bool SubsetCellFilter::lessThan(const CellFilterStub* other) const
{
  const SubsetCellFilter* S 
    = dynamic_cast<const SubsetCellFilter*>(other);

  TEST_FOR_EXCEPTION(S==0,
                     InternalError,
                     "argument " << other->toXML() 
                     << " to SubsetCellFilter::lessThan() should be "
                     "a SubsetCellFilter pointer.");

  return OrderedPair<CellFilter, CellPredicate>(superset_, predicate_)
    < OrderedPair<CellFilter, CellPredicate>(S->superset_, S->predicate_);
}

CellSet SubsetCellFilter::internalGetCells(const Mesh& mesh) const
{
  SUNDANCE_OUT(this->verb() > 1,
                   "SubsetCellFilter::internalGetCells()");
  CellSet super = superset_.getCells(mesh);

  int dim = superset_.dimension(mesh);

  CellType cellType = mesh.cellType(dim);

  predicate_.setMesh(mesh, dim);

  ExplicitCellSet* rtn = new ExplicitCellSet(mesh, dim, cellType);

  Set<int>& cells = rtn->cells();

  const CellPredicateBase* pred = predicate_.ptr().get();


  Array<int> cellLID;

  cellLID.reserve(mesh.numCells(dim));

  for (CellIterator i=super.begin(); i != super.end(); i++)
    {
      cellLID.append(*i);
    }

  Array<int> testResults(cellLID.size());
  pred->testBatch(cellLID, testResults);

  for (int i=0; i<cellLID.size(); i++)
    {
      SUNDANCE_OUT(this->verb() > 2,
                   "SubsetCellFilter is testing " << cellLID[i]);
      if (testResults[i]) 
        {
          SUNDANCE_OUT(this->verb() > 2,
                       "accepted " << cellLID[i]);
          cells.insert(cellLID[i]);
        }
      else
        {
          SUNDANCE_OUT(this->verb() > 2,
                       "rejected " << cellLID[i]);
        }
    }

  return rtn;
}
