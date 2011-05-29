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

#include "TSFPartitionedMatrixFactory.hpp"
#include "TSFVectorType.hpp"
#include "TSFLoadableBlockOperatorDecl.hpp"
#include "Teuchos_MPIComm.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorSpaceImpl.hpp"
#include "TSFSequentialIteratorImpl.hpp"
#include "TSFLoadableBlockOperatorImpl.hpp"
#endif

using namespace TSFExtended;
using namespace Teuchos;

PartitionedMatrixFactory::PartitionedMatrixFactory(
  const VectorSpace<double>& domain,
  int lowestLocalCol,
  const RCP<Array<int> >& isBCCol,
  const RCP<std::set<int> >& remoteBCCols,
  const VectorType<double>& domainVecType,
  const VectorSpace<double>& range,
  int lowestLocalRow,
  const RCP<Array<int> >& isBCRow,
  const VectorType<double>& rangeVecType
  )
  : 
  domain_(domain),
  internalDomain_(domain.getBlock(0)),
  bcDomain_(domain.getBlock(1)),
  isBCCol_(isBCCol),
  remoteBCCols_(remoteBCCols),
  domainVecType_(domainVecType),
  lowestLocalCol_(lowestLocalCol),
  highestLocalCol_(-1),
  range_(range),
  internalRange_(range.getBlock(0)),
  bcRange_(range.getBlock(1)),
  isBCRow_(isBCRow),
  rangeVecType_(rangeVecType),
  lowestLocalRow_(lowestLocalRow),
  highestLocalRow_(-1),
  blockFactory_(2),
  blockICMF_(2)
{
  highestLocalCol_ = lowestLocalCol_ + domain.numLocalElements();
  highestLocalRow_ = lowestLocalRow_ + range.numLocalElements();

  blockFactory_[0].resize(2);
  blockFactory_[1].resize(2);
  
  blockFactory_[0][0] = rangeVecType_.createMatrixFactory(internalDomain_, internalRange_);

  blockFactory_[0][1] = rangeVecType_.createMatrixFactory(bcDomain_, internalRange_);

  blockFactory_[1][1] = rangeVecType_.createMatrixFactory(bcDomain_, bcRange_);

  for (int i=0; i<2; i++)
  {
    blockICMF_[i].resize(2);
    for (int j=0; j<2; j++)
    {
      if (i==1 && j==0) continue;
      IncrementallyConfigurableMatrixFactory* icmf 
        = dynamic_cast<IncrementallyConfigurableMatrixFactory*>(blockFactory_[i][j].get());
      TEST_FOR_EXCEPTION(icmf==0, std::runtime_error,
        "block(" << i << ", " << j << ") is not an ICMF");
      blockICMF_[i][j] = icmf;
    }
  }
  
}

void PartitionedMatrixFactory::initializeNonzerosInRow(int globalRowIndex,
  int nElemsToInsert,
  const int* globalColumnIndices)
{
  if (globalRowIndex < lowestLocalRow_ || globalRowIndex >= highestLocalRow_) return;

  Array<int> bcCols;
  Array<int> intCols;
  bcCols.reserve(nElemsToInsert);
  intCols.reserve(nElemsToInsert);

  const Array<int>& isBCCol = *isBCCol_;

  for (int i=0; i<nElemsToInsert; i++)
  {
    int g = globalColumnIndices[i];
    if (g < lowestLocalCol_ || g >= highestLocalCol_)
    {
      if (remoteBCCols_->find(g) != remoteBCCols_->end()) bcCols.append(g);
      else intCols.append(g);
    }
    else
    {
      int localCol = g - lowestLocalCol_;
      if (isBCCol[localCol]) bcCols.append(g);
      else intCols.append(g);
    }
  }


  if ((*isBCRow_)[globalRowIndex - lowestLocalRow_])
  {
    if (intCols.size() > 0) /* do (BC, internal) block */
    {
      TEST_FOR_EXCEPTION(true, std::logic_error,
        "There should be no entries in the (BC, internal) block");
    }
    if (bcCols.size() > 0) /* do (BC, BC) block */
    {
      blockICMF_[1][1]->initializeNonzerosInRow(globalRowIndex, 
        bcCols.size(), &(bcCols[0]));
    }
  }
  else
  {
    if (intCols.size() > 0) /* do (internal, internal) block */
    {
      blockICMF_[0][0]->initializeNonzerosInRow(globalRowIndex, 
        intCols.size(), &(intCols[0]));
    }
    if (bcCols.size() > 0) /* do (internal, BC) block */
    {
      blockICMF_[0][1]->initializeNonzerosInRow(globalRowIndex, 
        bcCols.size(), &(bcCols[0]));
    }
  }
}


void PartitionedMatrixFactory::finalize() 
{
  blockICMF_[0][0]->finalize();
  blockICMF_[1][1]->finalize();
  blockICMF_[0][1]->finalize();
}


LinearOperator<double> PartitionedMatrixFactory::createMatrix() const
{

  RCP<LinearOpBase<double> > op 
    = rcp(new TSFExtended::LoadableBlockOperator<double>(domain_, lowestLocalCol_, isBCCol_, remoteBCCols_, range_, lowestLocalRow_, isBCRow_));
  LinearOperator<double> A = op;

  LinearOperator<double> A_ii = blockFactory_[0][0]->createMatrix();
  LinearOperator<double> A_ib = blockFactory_[0][1]->createMatrix();
  LinearOperator<double> A_bb = blockFactory_[1][1]->createMatrix();

  A.setBlock(0,0,A_ii);
  A.setBlock(0,1,A_ib);
  A.setBlock(1,1,A_bb);

  A.endBlockFill();

  return A;
}
