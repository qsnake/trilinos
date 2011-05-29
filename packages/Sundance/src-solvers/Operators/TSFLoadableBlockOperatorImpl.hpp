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

#ifndef TSFLOADABLEBLOCKOPERATOR_IMPL_HPP
#define TSFLOADABLEBLOCKOPERATOR_IMPL_HPP

#include "SundanceDefs.hpp"
#include "TSFLoadableBlockOperatorDecl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFSimpleBlockOpImpl.hpp"
#endif



namespace TSFExtended
{
using namespace Teuchos;

template <class Scalar> inline
LoadableBlockOperator<Scalar>:: LoadableBlockOperator(
  const VectorSpace<Scalar>& domain,
  int lowestLocalCol,
  const RCP<Array<int> >& isBCCol,
  const RCP<std::set<int> >& remoteBCCols,
  const VectorSpace<Scalar>& range,
  int lowestLocalRow,
  const RCP<Array<int> >& isBCRow)
  : SimpleBlockOp<Scalar>(domain, range),
    isBCCol_(isBCCol),
    isBCRow_(isBCRow),
    remoteBCCols_(remoteBCCols),
    lowestLocalRow_(lowestLocalRow),
    lowestLocalCol_(lowestLocalCol),
    highestLocalRow_(lowestLocalRow + range.numLocalElements()),
    highestLocalCol_(lowestLocalCol + domain.numLocalElements())
{
}
      
template <class Scalar> inline
void LoadableBlockOperator<Scalar>::addToRow(int globalRowIndex,
  int nElemsToInsert,
  const int* globalColumnIndices,
  const Scalar* elementValues) 
{
  if (globalRowIndex < lowestLocalRow_ || globalRowIndex >= highestLocalRow_) return;
  Array<int> bcCols;
  Array<int> intCols;
  Array<Scalar> bcVals;
  Array<Scalar> intVals;
  bcCols.reserve(nElemsToInsert);
  intCols.reserve(nElemsToInsert);
  bcVals.reserve(nElemsToInsert);
  intVals.reserve(nElemsToInsert);
        
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
      if (isBCCol[localCol]) 
      {
        bcCols.append(g);
        bcVals.append(elementValues[i]);
      }
      else 
      {
        intCols.append(g);
        intVals.append(elementValues[i]);
      }
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
      loadableBlock(1,1)->addToRow(globalRowIndex, bcCols.size(), 
        &(bcCols[0]), &(bcVals[0]));
    }
  }
  else
  {
    if (intCols.size() > 0) /* do (internal, internal) block */
    {
      loadableBlock(0,0)->addToRow(globalRowIndex, intCols.size(), 
        &(intCols[0]), &(intVals[0]));
    }
    if (bcCols.size() > 0) /* do (internal, BC) block */
    {
      loadableBlock(0,1)->addToRow(globalRowIndex, bcCols.size(), 
        &(bcCols[0]), &(bcVals[0]));
    }
  }
}


template <class Scalar> inline
void  LoadableBlockOperator<Scalar>::zero() 
{
  for (int i=0; i<this->numBlockRows(); i++)
  {
    for (int j=0; j<this->numBlockCols(); j++)
    {
      if (i==1 && j==0) continue;
      this->loadableBlock(i,j)->zero();
    }
  }
}

    

template <class Scalar> inline
RCP<LoadableMatrix<Scalar> > 
LoadableBlockOperator<Scalar>::loadableBlock(int i, int j)
{
  return this->getNonconstBlock(i,j).matrix();
}


}

#endif
