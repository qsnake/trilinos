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

#include "TSFEpetraMatrixFactory.hpp"
#include "TSFEpetraMatrix.hpp"
#include "TSFEpetraVector.hpp"
#include "TSFVectorSpaceDecl.hpp"  
#include "TSFVectorDecl.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_MPIComm.hpp"
#include "TSFLinearOperatorDecl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#include "TSFVectorImpl.hpp"
#endif


using namespace TSFExtended;
using namespace Teuchos;
using namespace Thyra;

EpetraMatrixFactory::EpetraMatrixFactory(const RCP<const EpetraVectorSpace>& domain,
                         const RCP<const EpetraVectorSpace>& range)
  : graph_(rcp(new Epetra_CrsGraph(Copy, *(range->epetraMap()), 0))),
    range_(range),
    domain_(domain)
{}


void EpetraMatrixFactory::finalize()
{
  int ierr = graph_->FillComplete(*(domain_->epetraMap()), *(range_->epetraMap()));

  TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error, 
                     "EpetraMatrixFactory::finalize() failed during call "
                     "to FillComplete(). Error code was " << ierr);

  if (!graph_->StorageOptimized())
    {
      ierr = graph_->OptimizeStorage();
      
      TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error, 
                         "EpetraMatrixFactory::freezeValues() failed during call "
                         "to OptimizeStorage(). Error code was " << ierr);
    }
}

void EpetraMatrixFactory::initializeNonzerosInRow(int globalRowIndex,
                                          int nElemsToInsert,
                                          const int* globalColumnIndices)
{
  int ierr = graph_->InsertGlobalIndices(globalRowIndex,
                                         nElemsToInsert,
                                         (int*) globalColumnIndices);
  
  TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error, 
                     "failed to add to row " << globalRowIndex
                     << " in EpetraMatrixFactory::setRowValues() with nnz="
                     << nElemsToInsert 
                     << ". Error code was " << ierr);
}


void EpetraMatrixFactory::initializeNonzeroBatch(int numRows, 
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
          graph_->InsertGlobalIndices(globalRowIndices[row], 
                                      numColumnsPerRow, (int*) cols);
        }
    }
}


void EpetraMatrixFactory::configure(int lowestRow,
                            const std::vector<int>& rowPtrs,
                            const std::vector<int>& nnzPerRow,
                            const std::vector<int>& data)
{

  graph_ = rcp(new Epetra_CrsGraph(Copy, *(range_->epetraMap()),
                                   (const int*) &(nnzPerRow[0]),
                                   true));
  
  for (unsigned int i=0; i<rowPtrs.size(); i++)
    {
      graph_->InsertGlobalIndices(lowestRow + i, nnzPerRow[i],
                                  (int*) &(data[rowPtrs[i]]));
    }

  finalize();
}

const Epetra_CrsGraph& EpetraMatrixFactory::graph() const 
{
  return *(graph_.get());
}


LinearOperator<double> EpetraMatrixFactory::createMatrix() const
{
  RCP<LinearOpBase<double> > A 
    = rcp(new EpetraMatrix(graph(), epDomain(), epRange()));
  return A;
}
