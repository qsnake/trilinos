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
#include "TSFEpetraMatrixMatrixProduct.hpp"
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
using Sundance::RuntimeError;


LinearOperator<double> epetraLeftScale(
  const Vector<double>& d,
  const LinearOperator<double>& A)
{
  /* Extract the underlying Epetra matrix. Type checking is done
   * ny rcp_dynamic_cast, so we need no error check here. */
  RCP<const Epetra_CrsMatrix> A_crs = EpetraMatrix::getConcretePtr(A);
  
  /* Make a deep copy of A */
  RCP<Epetra_CrsMatrix> mtxCopy = rcp(new Epetra_CrsMatrix(*A_crs));

  /* Extract the underlying Epetra vector. Type checking is done
   * internally, so we need no error check here. */
  const Epetra_Vector& epv = EpetraVector::getConcrete(d);
  
  /* Scale the copy */
  mtxCopy->LeftScale(epv);

  /* Prepare an operator object for the scaled matrix */
  RCP<const EpetraVectorSpace> domain 
    = rcp_dynamic_cast<const EpetraVectorSpace>(A.domain().ptr());

  RCP<const EpetraVectorSpace> range 
    = rcp_dynamic_cast<const EpetraVectorSpace>(A.range().ptr());

  RCP<LinearOpBase<double> > rtn 
    = rcp(new EpetraMatrix(mtxCopy, domain, range));
  return rtn;
  
}

LinearOperator<double> epetraRightScale(
  const LinearOperator<double>& A,
  const Vector<double>& d)
{
  /* Extract the underlying Epetra matrix. Type checking is done
   * ny rcp_dynamic_cast, so we need no error check here. */
  RCP<const Epetra_CrsMatrix> A_crs = EpetraMatrix::getConcretePtr(A);
  
  /* Make a deep copy of A */
  RCP<Epetra_CrsMatrix> mtxCopy = rcp(new Epetra_CrsMatrix(*A_crs));

  /* Extract the underlying Epetra vector. Type checking is done
   * internally, so we need no error check here. */
  const Epetra_Vector& epv = EpetraVector::getConcrete(d);
  
  /* Scale the copy */
  mtxCopy->RightScale(epv);

  /* Prepare an operator object for the scaled matrix */
  RCP<const EpetraVectorSpace> domain 
    = rcp_dynamic_cast<const EpetraVectorSpace>(A.domain().ptr());

  RCP<const EpetraVectorSpace> range 
    = rcp_dynamic_cast<const EpetraVectorSpace>(A.range().ptr());

  RCP<LinearOpBase<double> > rtn 
    = rcp(new EpetraMatrix(mtxCopy, domain, range));
  return rtn;
  
}


LinearOperator<double> epetraMatrixMatrixProduct(
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
  

  /* Get the row map from A. We will need this to build the target matrix C */
  const Epetra_Map* rowmap 
    = transA ? &(A_crs->DomainMap()) : &(A_crs->RowMap());

  /* make the target matrix */
  RCP<Epetra_CrsMatrix> C = rcp(new Epetra_CrsMatrix(Copy, *rowmap, 1));

  /* Carry out the multiplication */
  int ierr 
    = EpetraExt::MatrixMatrix::Multiply(*A_crs, transA, *B_crs, transB, *C);
  TEST_FOR_EXCEPTION(ierr != 0, RuntimeError,
    "EpetraExt Matrix-matrix multiply failed with error code ierr=" << ierr);

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
