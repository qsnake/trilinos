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

#ifndef SUNDANCE_MATRIXVECTORASSEMBLYKERNEL_H
#define SUNDANCE_MATRIXVECTORASSEMBLYKERNEL_H

#include "SundanceDefs.hpp"
#include "SundanceVectorFillingAssemblyKernel.hpp"

namespace Sundance
{
using namespace Teuchos;

/** 
 * MatrixVectorAssemblyKernel does assembly of a matrix and vector
 */
class MatrixVectorAssemblyKernel : public VectorFillingAssemblyKernel
{
public:
  /** */
  MatrixVectorAssemblyKernel(
    const Array<RCP<DOFMapBase> >& rowMap,
    const Array<RCP<Array<int> > >& isBCRow,
    const Array<int>& lowestLocalRow,
    const Array<RCP<DOFMapBase> >& colMap,
    const Array<RCP<Array<int> > >& isBCCol,
    const Array<int>& lowestLocalCol,
    LinearOperator<double> A,
    Array<Vector<double> > b,
    bool partitionBCs,
    int verb)
    : VectorFillingAssemblyKernel(rowMap, isBCRow, lowestLocalRow, 
      b, partitionBCs, verb),
      mat_(rowMap.size()),
      cmb_(colMap, isBCCol, lowestLocalCol, partitionBCs, verb)
    {
      init(rowMap, colMap, A, partitionBCs);
    }

  /** */
  void prepareForWorkSet(
    const Array<Set<int> >& requiredTests,
    const Array<Set<int> >& requiredUnks,
    RCP<StdFwkEvalMediator> mediator) ;

  /** */
  void fill(bool isBC,
    const IntegralGroup& group,
    const RCP<Array<double> >& localValues) ;

protected:

  /** */
  void init(
  const Array<RCP<DOFMapBase> >& rowMap,
  const Array<RCP<DOFMapBase> >& colMap,
  LinearOperator<double> A,
  bool partitionBCs);

  /** */
  void writeLSMs(int blockRow, int blockCol,
    bool useCofacetCells,
    int numTestNodes, 
    int nTestFuncs, 
    int testFuncIndex, 
    const Array<int>& rowDof,
    int numUnkNodes, 
    int nUnkFuncs, 
    int unkFuncIndex, 
    const Array<int>& colDof,
    const Array<double>& localValues) const ;

  /** */
  void insertLocalMatrixBatch(
    bool isBCRqc,
    bool useCofacetCells,
    const Array<int>& testID, 
    const Array<int>& testBlock, 
    const Array<int>& unkID,
    const Array<int>& unkBlock,
    const Array<double>& localValues) const ;

protected:
  const MapBundle& rmb() const {return mapBundle();}
  const MapBundle& cmb() const {return cmb_;}

private:
  LinearOperator<double> A_;
  Array<Array<LoadableMatrix<double>* > > mat_;
  mutable MapBundle cmb_;
};

}



#endif
