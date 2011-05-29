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

#ifndef TSFPARTITIONEDMATRIXFACTORYDECL_HPP
#define TSFPARTITIONEDMATRIXFACTORYDECL_HPP

#include "TSFIncrementallyConfigurableMatrixFactory.hpp"
#include "TSFCollectivelyConfigurableMatrixFactory.hpp"
#include "TSFMatrixFactory.hpp"
#include "TSFVectorSpaceDecl.hpp"
#include "TSFVectorType.hpp"
#include "SundanceHandleable.hpp"
#include "SundancePrintable.hpp"
#include "Epetra_CrsGraph.h"
#include <set>

namespace TSFExtended
{
  using namespace Teuchos;
  using namespace Thyra;

  /** */
  class PartitionedMatrixFactory 
    : public MatrixFactory<double>,
      public IncrementallyConfigurableMatrixFactory
  {
  public:

    /** Construct an uninitialized PartitionedMatrixFactory */
    PartitionedMatrixFactory(
      const VectorSpace<double>& domain,
      int lowestLocalCol,
      const RCP<Array<int> >& isBCCol,
      const RCP<std::set<int> >& remoteBCCols,
      const VectorType<double>& domainVecType,
      const VectorSpace<double>& range,
      int lowestLocalRow,
      const RCP<Array<int> >& isBCRow,
      const VectorType<double>& rangeType
      );

    /** */
    virtual ~PartitionedMatrixFactory(){;}

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

    /** */
    void finalize();

    /** */
    virtual LinearOperator<double> createMatrix() const ;

  protected:

  private:
    VectorSpace<double> domain_;
    VectorSpace<double> internalDomain_;
    VectorSpace<double> bcDomain_;
    RCP<Array<int> > isBCCol_;
    RCP<std::set<int> > remoteBCCols_;
    VectorType<double> domainVecType_;
    int lowestLocalCol_;
    int highestLocalCol_;
    VectorSpace<double> range_;
    VectorSpace<double> internalRange_;
    VectorSpace<double> bcRange_;
    RCP<Array<int> > isBCRow_;
    VectorType<double> rangeVecType_;
    int lowestLocalRow_;
    int highestLocalRow_;


    Array<Array<RCP<MatrixFactory<double> > > > blockFactory_;
    Array<Array<IncrementallyConfigurableMatrixFactory*> > blockICMF_;
  };
}


#endif
