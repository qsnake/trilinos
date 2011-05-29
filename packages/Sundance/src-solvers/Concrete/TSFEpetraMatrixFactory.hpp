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

#ifndef TSFEPETRAMATRIXFACTORY_HPP
#define TSFEPETRAMATRIXFACTORY_HPP

#include "TSFEpetraVectorSpace.hpp"
#include "TSFIncrementallyConfigurableMatrixFactory.hpp"
#include "TSFCollectivelyConfigurableMatrixFactory.hpp"
#include "TSFMatrixFactory.hpp"
#include "SundanceHandleable.hpp"
#include "SundancePrintable.hpp"
#include "Epetra_CrsGraph.h"

namespace TSFExtended
{
  using namespace Teuchos;
  using namespace Thyra;

  /** */
  class EpetraMatrixFactory : public MatrixFactory<double>,
                      public IncrementallyConfigurableMatrixFactory,
                      public CollectivelyConfigurableMatrixFactory
  {
  public:

    /** Construct an uninitialized EpetraMatrixFactory */
    EpetraMatrixFactory(const RCP<const EpetraVectorSpace>& domain,
                const RCP<const EpetraVectorSpace>& range);

    /** */
    const RCP<const EpetraVectorSpace>& epRange() const {return range_;}

    /** */
    const RCP<const EpetraVectorSpace>& epDomain() const {return domain_;}


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
                                         const int* globalColumnIndices) ;

    /** 
     * Initialize nonzeros in a batch of rows. 
     */
    virtual void initializeNonzeroBatch(int numRows, 
                                        int rowBlockSize,
                                        const int* globalRowIndices,
                                        int numColumnsPerRow,
                                        const int* globalColumnIndices,
                                        const int* skipRow);

    /** Configure all rows at once */
    virtual void configure(int lowestRow,
                           const std::vector<int>& rowPtrs,
                           const std::vector<int>& nnzPerRow,
                           const std::vector<int>& data);

    /** */
    void finalize();

    /** */
    const Epetra_CrsGraph& graph() const ;

    /** */
    virtual LinearOperator<double> createMatrix() const ;

  protected:

  private:

    /** */
    RCP<Epetra_CrsGraph> graph_;

    /** */
    RCP<const EpetraVectorSpace> range_;

    /** */
    RCP<const EpetraVectorSpace> domain_;
  };
}

#endif
