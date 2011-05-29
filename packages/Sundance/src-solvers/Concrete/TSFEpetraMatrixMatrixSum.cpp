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
#include "TSFEpetraMatrixMatrixSum.hpp"
#include "TSFEpetraVector.hpp"
#include "SundanceExceptions.hpp"
#include "EpetraExt_MatrixMatrix.h"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorImpl.hpp"
#include "TSFLinearOperatorImpl.hpp"
#endif



namespace TSFExtended
{
using namespace Teuchos;


LinearOperator<double> epetraMatrixMatrixSum(
  const LinearOperator<double>& A,
  const LinearOperator<double>& B)
{
  /* Extract the underlying Epetra matrix for A. Type checking is done
   * ny rcp_dynamic_cast, so we need no error check here. */
  RCP<const Epetra_CrsMatrix> A_crs = EpetraMatrix::getConcretePtr(A);

  /* Extract the underlying Epetra matrix for A. Type checking is done
   * ny rcp_dynamic_cast, so we need no error check here. */
  RCP<const Epetra_CrsMatrix> B_crs = EpetraMatrix::getConcretePtr(B);
  
  bool transA = false;
  bool transB = false;

  TEST_FOR_EXCEPTION(A.range() != B.range(), RuntimeError,
    "incompatible operand ranges in epetraMatrixMatrixSum()"
    << std::endl << "A.range()=" << A.range()
    << std::endl << "B.range()=" << B.range()
    );
  

  TEST_FOR_EXCEPTION(A.domain() != B.domain(), RuntimeError,
    "incompatible operand domains in epetraMatrixMatrixSum()"
    << std::endl << "A.domain()=" << A.domain()
    << std::endl << "B.domain()=" << B.domain()
    );
  

  /* Get the row map from A. We will need this to build the target matrix C */
  const Epetra_Map* rowmap 
    = transA ? &(A_crs->DomainMap()) : &(A_crs->RowMap());

  /* make the target matrix */
  RCP<Epetra_CrsMatrix> C = rcp(new Epetra_CrsMatrix(Copy, *rowmap, 1));
  Epetra_CrsMatrix* CPtr = C.get();

  /* Carry out the multiplication */
  int ierr 
    = EpetraExt::MatrixMatrix::Add(
      *A_crs, transA, 1.0, 
      *B_crs, transB, 1.0, CPtr);
  TEST_FOR_EXCEPTION(ierr != 0, RuntimeError,
    "EpetraExt Matrix-matrix add failed with error code ierr=" << ierr);

  /* Need to call FillComplete on the result */
  C->FillComplete();

  /* Prepare an operator object for the scaled matrix */
  RCP<const EpetraVectorSpace> range 
    = rcp_dynamic_cast<const EpetraVectorSpace>(A.range().ptr());

  RCP<const EpetraVectorSpace> domain 
    = rcp_dynamic_cast<const EpetraVectorSpace>(B.domain().ptr());

  RCP<LinearOpBase<double> > rtn 
    = rcp(new EpetraMatrix(C, domain, range));
  return rtn;
  
}

}
