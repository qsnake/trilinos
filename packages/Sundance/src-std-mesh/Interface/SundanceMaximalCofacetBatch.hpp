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

#ifndef SUNDANCE_MAXIMALCOFACETBATCH_H
#define SUNDANCE_MAXIMALCOFACETBATCH_H

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace Sundance {
using namespace Teuchos;

/** 
 * MaximalCofacetBatch is used to store the maximal cofacets
 * of a batch of cells of codimension one. Interior cells will
 * have two cofacets, while boundary cells will have only one. Every 
 * cell in the batch required to have the same number of cofacets, which
 * can be arranged at the user level by careful distinction between
 * external boundaries and internal boundaries. 
 */
class MaximalCofacetBatch
{
public:
  /** 
   * Initialize an empty batch
   */
  MaximalCofacetBatch();

  /** 
   * Change the number of cells in the batch.
   */
  void reset(int numCells);

  /** 
   * Change the number of cells and cofacets in the batch.
   */
  void reset(int numCells, int numCofacets);

  /** 
   * Return the number of cells in the batch
   */
  int numCells() const {return numCells_;}

  /** 
   * 
   */
  int numCofacets() const {return numCofacets_;}

  /** 
   * Return the LID of the n-th cofacet of the c-th cell in the batch. 
   * 
   * \param c The index (within the batch) of the cell whose cofacets are
   * requested.
   * \param n The number of the cofacet requested. This can only be 0 or 1,
   * and must be less than numCells().
   * \param facetIndex The c-th cell in the batch is one of the facets of
   * its maximal cofacets. Its facet index within that cell is returned
   * via reference as facetIndex.
   * \returns The LID of the requested cofacet.
   */
  int cofacetLID(int c, int n, int& facetIndex) const ; 


  /** 
   * Pick one specified cofacet for each cell in the batch.
   *
   * \param cofacetNumbers the c-th entry selects cofacet for cell c. Each
   * entry in the array should be either 0 or 1.
   * \param cofacets array for storage of the results.
   * \param facet index array for storage of the results.
   */
  void getSpecifiedCofacets(const Array<int>& cofacetNumbers,
    RCP<Array<int> >& cofacets, RCP<Array<int> >& facetIndices) const ;


  /** 
   * Pick one specified cofacet for each cell in the batch.
   *
   * \param cofacetNumber which cofacet to select
   * \param cofacets array for storage of the results.
   * \param facet index array for storage of the results.
   */
  void getSpecifiedCofacets(int cofacetNumber,
    RCP<Array<int> >& cofacets, RCP<Array<int> >& facetIndices) const ;

  /** 
   * Add a cell with a single cofacet
   */
  void addSingleCofacet(int c, int cofacetLID, int facetIndex) ;

  /** 
   * Add a cell with two cofacets
   */
  void addTwoCofacets(int c, 
    int cofacet1, int facetIndex1,
    int cofacet2, int facetIndex2);

private:
  Array<RCP<Array<int> > > cofacetLIDs_;

  Array<RCP<Array<int> > > facetIndices_;

  int numCells_;

  int numCofacets_;
};

}

#endif


