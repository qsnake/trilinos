/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// **********************************************************************/
/* @HEADER@ */

#ifndef TSF_INCREMENTALLYCONFIGURABLEMATRIXFACTORY_HPP
#define TSF_INCREMENTALLYCONFIGURABLEMATRIXFACTORY_HPP

#include "SundanceDefs.hpp"

namespace TSFExtended
{
  /** 
   * Class IncrementallyConfigurableMatrixFactory provides an abstract 
   * interface for row-at-a-time configuration of matrix factories.
   */
  class IncrementallyConfigurableMatrixFactory
  {
  public:
    /** Virtual dtor */
    virtual ~IncrementallyConfigurableMatrixFactory(){;}

    /** Initialize a set of nonzero elements in the matrix's graph.
     * @param globalRowIndex the global index of the row to which these
     * elements belong.
     * @param nElemsToInsert the number of elements being inserted in this
     * step
     * @param globalColumnIndices array of column indices. Must 
     * be nElemsToInsert in length. 
     */
    virtual void initializeNonzerosInRow(int globalRowIndex,
                                         int nElemsToInsert,
                                         const int* globalColumnIndices) = 0 ;

    /** 
     * Initialize nonzeros in a batch of rows. 
     */
    virtual void initializeNonzeroBatch(int numRows, 
                                        int rowBlockSize,
                                        const int* globalRowIndices,
                                        int numColumnsPerRow,
                                        const int* globalColumnIndices,
                                        const int* skipRow);

    /** Finalize values of the matrix. This is a hook for any
     * implementation-dependent steps that must be done after
     * loading of elements. */
    virtual void finalize() = 0 ;

  private:
    
    
  };

  /* Default implementation of initializeElementBatch */
  inline void IncrementallyConfigurableMatrixFactory
  ::initializeNonzeroBatch(int numRows, 
                           int rowBlockSize,
                           const int* globalRowIndices,
                           int numColumnsPerRow,
                           const int* globalColumnIndices,
                           const int* skipRow)
  {
    int numRowBlocks = numRows/rowBlockSize;
    int row = 0;

    for (int rb=0; rb<numRowBlocks; rb++)
      {
        const int* cols = globalColumnIndices + rb*numColumnsPerRow;
        for (int r=0; r<rowBlockSize; r++, row++)
          {
            if (skipRow[row]) continue;
            initializeNonzerosInRow(globalRowIndices[row], 
                                    numColumnsPerRow, cols);
          }
      }
  }
}

#endif
