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

#include "TSFPoissonBoltzmannJacobian.hpp"
#include "TSFEpetraMatrix.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorImpl.hpp"
#include "TSFLinearOperatorImpl.hpp"
#endif
using namespace TSFExtended;
using namespace Teuchos;


PoissonBoltzmannJacobian
::PoissonBoltzmannJacobian(int nLocalRows, 
                           const VectorType<double>& type)
  : OperatorBuilder<double>(nLocalRows, type), op_(), nLocalRows_(nLocalRows),
    h_(1.0)
{
  h_ = 1.0/((double) domain().dim() - 1);
}

void PoissonBoltzmannJacobian::setEvalPoint(const Vector<double>& x)
{
  
  int rank = MPIComm::world().getRank();
  int nProc = MPIComm::world().getNProc();
  RCP<MatrixFactory<double> > mFact 
    = vecType().createMatrixFactory(domain(), range());
  
  int lowestLocalRow = nLocalRows_ * rank;

  IncrementallyConfigurableMatrixFactory* icmf 
    = dynamic_cast<IncrementallyConfigurableMatrixFactory*>(mFact.get());
  for (int i=0; i<nLocalRows_; i++)
    {
      int row = lowestLocalRow + i;
      Array<int> colIndices;
      if ((rank==0 && i==0) || (rank==nProc-1 && i==nLocalRows_-1))
        {
          colIndices = tuple(row);
        }
      else
        {
          colIndices = tuple(row-1, row, row+1);
        }
      icmf->initializeNonzerosInRow(row, colIndices.size(),
                                    &(colIndices[0]));
    }
  icmf->finalize();
      
  op_ = mFact->createMatrix();
      
  RCP<LoadableMatrix<double> > mat = op_.matrix();

  /* fill in with the Laplacian operator plus exp(-x) */
  for (int i=0; i<nLocalRows_; i++)
    {
      int row = lowestLocalRow + i;
      Array<int> colIndices;
      Array<double> colVals;
      if ((rank==0 && i==0) || (rank==nProc-1 && i==nLocalRows_-1))
        {
          colIndices = tuple(row);
          colVals = tuple(1.0);
        }
      else
        {
          colIndices = tuple(row-1, row, row+1);
          colVals = tuple(1.0/h_/h_, 
                          -2.0/h_/h_ + exp(-x.getElement(row)), 
                          1.0/h_/h_);
        }
      mat->addToRow(row, colIndices.size(), 
                    &(colIndices[0]), &(colVals[0]));
    }
}
