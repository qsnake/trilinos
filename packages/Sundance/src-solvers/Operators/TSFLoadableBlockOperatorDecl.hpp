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

#ifndef TSFLOADABLEBLOCKOPERATOR_DECL_HPP
#define TSFLOADABLEBLOCKOPERATOR_DECL_HPP

#include "SundanceDefs.hpp"
#include "Teuchos_RCP.hpp"
#include "TSFSimpleBlockOpDecl.hpp"
#include "TSFLoadableMatrix.hpp"

namespace TSFExtended
{
using namespace Teuchos;
  /** 
   * Class LoadableBlockOperator provides a LoadableMatrix interface
   * for a physically-partitioned block 2x2 matrix, making it appear
   * to the fill routine as if the block matrix is a single matrix. This 
   * is intended for filling systems where the internal and BC equations
   * and unknowns are stored in physically separate blocks.
   */
  template <class Scalar>
  class LoadableBlockOperator 
    : public SimpleBlockOp<Scalar>,
      public LoadableMatrix<Scalar>
  {
  public:
    /** */
    LoadableBlockOperator(
      const VectorSpace<Scalar>& domain,
      int lowestLocalCol,
      const RCP<Array<int> >& isBCCol,
      const RCP<std::set<int> >& remoteBCCols,
      const VectorSpace<Scalar>& range,
      int lowestLocalRow,
      const RCP<Array<int> >& isBCRow);
      
    /** Virtual dtor */
    virtual ~LoadableBlockOperator(){;}

    /** Insert a set of elements in a row, adding to any previously
     * existing values.  The nonzero structure of the matrix must have
     * been determined at construction time. 
     *
     * @param globalRowIndex the global index of the row to which these
     * elements belong.
     * @param nElemsToInsert the number of elements being inserted in this
     * step
     * @param globalColumnIndices array of column indices. Must 
     * be nElemsToInsert in length. 
     * @param elements array of element values. Must be nElemsToInsert in
     * length
     */
    virtual void addToRow(int globalRowIndex,
      int nElemsToInsert,
      const int* globalColumnIndices,
      const Scalar* elementValues) ;

    /** Set all elements to zero, preserving the existing structure */
    virtual void zero() ;

    
  private:
    

    /** Cast a block to LoadableMatrix, with safety checks */
    RCP<LoadableMatrix<Scalar> > loadableBlock(int i, int j);

    RCP<Array<int> > isBCCol_;
    RCP<Array<int> > isBCRow_;
    RCP<std::set<int> > remoteBCCols_;
    int lowestLocalRow_;
    int lowestLocalCol_;
    int highestLocalRow_;
    int highestLocalCol_;
  };
}

#endif
