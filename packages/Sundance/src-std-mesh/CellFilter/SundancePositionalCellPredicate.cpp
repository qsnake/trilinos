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

#include "SundancePositionalCellPredicate.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

bool PositionalCellPredicate::lessThan(const CellPredicateBase* other) const
{
  TEST_FOR_EXCEPTION(dynamic_cast<const PositionalCellPredicate*>(other) == 0,
                     InternalError,
                     "argument " << other->toXML() 
                     << " to PositionalCellPredicate::lessThan() should be "
                     "a PositionalCellPredicate pointer.");

  return func_.get() < dynamic_cast<const PositionalCellPredicate*>(other)->func_.get();
}

void PositionalCellPredicate::testBatch(const Array<int>& cellLID,
                                        Array<int>& results) const
{
  results.resize(cellLID.size());

  if (cellDim()==0)
    {
      for (int i=0; i<cellLID.size(); i++)
        {
          results[i] = (*func_)(mesh().nodePosition(cellLID[i]));
        }
    }
  else
    {
      Array<int> facetLIDs;
      Array<int> facetSigns;
      int nf = mesh().numFacets(cellDim(), cellLID[0], 0);
      mesh().getFacetLIDs(cellDim(), cellLID, 0, facetLIDs, facetSigns);
      for (int c=0; c<cellLID.size(); c++)
        {
          results[c] = true;
          for (int f=0; f<nf; f++)
            {
              int fLID = facetLIDs[c*nf + f];
              if ((*func_)(mesh().nodePosition(fLID)) == false)
                {
                  results[c] = false;
                  continue;
                }
            }
        }
    }
}

XMLObject PositionalCellPredicate::toXML() const 
{
  XMLObject rtn("PositionalCellPredicate");
  return rtn;
}

bool PointCellPredicateFunctor::operator()(const Point& x) const
{
  for (int i=0; i<x_.dim(); i++)
  {
    if (fabs(x[i]-x_[i]) > tol_) return false;
  }
  return true;
}


bool CoordinateValueCellPredicateFunctor::operator()(const Point& x) const
{
  return (fabs(x[direction_]-value_) < tol_);
}
