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

#include "SundanceMaximalCofacetBatch.hpp"
#include "SundanceExceptions.hpp"

using namespace Sundance;
using namespace Teuchos;
using namespace Sundance;

MaximalCofacetBatch::MaximalCofacetBatch()
  : cofacetLIDs_(2),
    facetIndices_(2),
    numCells_(-1),
    numCofacets_(-1)
{
  for (int i=0; i<2; i++)
  {
    cofacetLIDs_[i] = rcp(new Array<int>());
    facetIndices_[i] = rcp(new Array<int>());
  }
}

void MaximalCofacetBatch::reset(int numCells, int numCofacets)
{
  TEST_FOR_EXCEPTION(numCofacets < 0 || numCofacets > 2, RuntimeError,
    "invalid number of maximal cofacets = " << numCofacets);
  numCofacets_ = numCofacets;
  reset(numCells);
}

void MaximalCofacetBatch::reset(int numCells)
{
  numCells_ = numCells;
  TEST_FOR_EXCEPTION(numCells_ <= 0, RuntimeError,
    "invalid number of cells = " << numCells_);
  for (int i=0; i<numCofacets_; i++)
  {
    cofacetLIDs_[i]->resize(numCells_);
    facetIndices_[i]->resize(numCells_);
  }
}

void MaximalCofacetBatch::addSingleCofacet(int c, 
  int cofacetLID, int facetIndex)
{
  TEST_FOR_EXCEPTION(numCofacets_ != 1, RuntimeError,
    "addSingleCofacet called for a batch configured with " << numCofacets_ 
    << " cofacets");

  cofacetLIDs_[0]->operator[](c) = cofacetLID;
  facetIndices_[0]->operator[](c) = facetIndex;
}
void MaximalCofacetBatch::addTwoCofacets(int c, 
  int cofacet1, int facetIndex1,
  int cofacet2, int facetIndex2
)
{
  TEST_FOR_EXCEPTION(numCofacets_ != 2, RuntimeError,
    "addTwoCofacets called for a batch configured with " << numCofacets_ 
    << " cofacets");

  cofacetLIDs_[0]->operator[](c) = cofacet1;
  facetIndices_[0]->operator[](c) = facetIndex1;

  cofacetLIDs_[1]->operator[](c) = cofacet2;
  facetIndices_[1]->operator[](c) = facetIndex2;
}

int MaximalCofacetBatch::cofacetLID(int c, int n, int& facetIndex) const
{
  TEST_FOR_EXCEPTION(n >= numCofacets_, RuntimeError,
    "invalid cofacet number n=" << n);
  TEST_FOR_EXCEPTION(c >= numCells_, RuntimeError,
    "invalid cell number c=" << c);

  facetIndex = facetIndices_[n]->operator[](c);
  return cofacetLIDs_[n]->operator[](c);
}

void MaximalCofacetBatch::getSpecifiedCofacets(
  const Array<int>& cofacetNumbers,
  RCP<Array<int> >& cofacets,
  RCP<Array<int> >& facetIndices) const
{
  TEST_FOR_EXCEPTION((int) cofacetNumbers.size() != numCells(),
    RuntimeError,
    "mismatch between cofacet batch size (" << numCells() << ") and "
    "requested number of cofacets (" << cofacetNumbers.size() << ")");

  cofacets->resize(numCells());
  facetIndices->resize(numCells());

  for (int c=0; c<numCells(); c++)
  {
    (*cofacets)[c] = cofacetLID(c, cofacetNumbers[c], (*facetIndices)[c]);
  }
}

void MaximalCofacetBatch::getSpecifiedCofacets(
  int cofacetNumber,
  RCP<Array<int> >& cofacets,
  RCP<Array<int> >& facetIndices) const
{
  TEST_FOR_EXCEPTION(cofacetNumber < 0 || cofacetNumber>1,    
    RuntimeError,
    "invalid cofacet number=" << cofacetNumber);

  cofacets = cofacetLIDs_[cofacetNumber];
  facetIndices = facetIndices_[cofacetNumber];
}
