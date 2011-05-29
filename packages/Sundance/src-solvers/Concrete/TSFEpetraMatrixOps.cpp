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

#include "TSFEpetraMatrix.hpp"
#include "TSFEpetraVector.hpp"
#include "TSFVectorSpaceDecl.hpp"  // changed from Impl
#include "TSFVectorDecl.hpp"
#include "TSFLinearOperatorDecl.hpp"  // changed from Impl

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorImpl.hpp"
#include "TSFLinearOperatorImpl.hpp"
#endif

using namespace TSFExtended;
using namespace Teuchos;
using namespace Thyra;



namespace TSFExtended
{

Vector<double> getEpetraDiagonal(const LinearOperator<double>& A)
{
  /* Extract the underlying Epetra matrix. Type checking is done
   * ny rcp_dynamic_cast, so we need no error check here. */
  RCP<const Epetra_CrsMatrix> A_crs = EpetraMatrix::getConcretePtr(A);

  VectorSpace<double> rowSpace = A.domain();
  Vector<double> rtn = rowSpace.createMember();

  Epetra_Vector* xPtr = EpetraVector::getConcretePtr(rtn);
  A_crs->ExtractDiagonalCopy(*xPtr);

  return rtn;
}


LinearOperator<double> makeEpetraDiagonalMatrix(const Vector<double>& d)
{
  VectorSpace<double> space = d.space();
  RCP<const EpetraVectorSpace> eps 
    = rcp_dynamic_cast<const EpetraVectorSpace>(space.ptr());

  EpetraMatrixFactory mf(eps, eps);

  int nLocal = space.numLocalElements();
  int offset = space.lowestLocallyOwnedIndex();
  for (int i=0; i<nLocal; i++)
  {
    int rowIndex = offset + i;
    mf.initializeNonzerosInRow(rowIndex, 1, &rowIndex);
  }

  mf.finalize();
  LinearOperator<double> rtn = mf.createMatrix();

  RCP<EpetraMatrix> epm = rcp_dynamic_cast<EpetraMatrix>(rtn.ptr());
  epm->zero();
  
  for (int i=0; i<nLocal; i++)
  {
    int rowIndex = offset + i;
    double val = d.getElement(rowIndex);
    epm->addToRow(rowIndex, 1, &rowIndex, &val);
  }
  
  return rtn;
}

}
