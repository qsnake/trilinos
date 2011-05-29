// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "DenseLinAlgLAPack.hpp"
#include "DenseLinAlgPack_LAPACK_Cpp.hpp"
#include "DenseLinAlgPack_DMatrixAsTriSym.hpp"
#include "Teuchos_TestForException.hpp"

namespace {
template< class T >
inline
T my_min( const T& v1, const T& v2 ) { return v1 < v2 ? v1 : v2; }
} // end namespace

void DenseLinAlgLAPack::potrf( DMatrixSliceTriEle* A )
{
  FortranTypes::f_int info;
  LAPACK_Cpp::potrf( A->uplo(), A->rows(), A->gms().col_ptr(1), A->gms().max_rows(), &info );
  if( info != 0 ) {
    TEST_FOR_EXCEPTION(
      info < 0, std::invalid_argument
      ,"potrf(...): Error, Invalid argument "
      << -info << " sent to LAPACK function xPOTRF(...)" );
    // info > 0
    TEST_FOR_EXCEPTION(
      true, FactorizationException
      ,"potrf(...): Error, Minor of order "
      << info << " is not positive definite, the factorization "
      "could not be completed" );
  }
  // If you get here all went okay and A has been factorized and now
  // hold the cholesky factor of input A.
}

void DenseLinAlgLAPack::geqrf( DMatrixSlice* A, DVectorSlice* tau, DVectorSlice* work )
{
  FortranTypes::f_int info;
  if( tau->dim() != my_min( A->rows(), A->cols() ) ) {
    TEST_FOR_EXCEPTION(
      true, std::invalid_argument, "geqrf(...): Error, tau is not sized correctly!" );
  }
  LAPACK_Cpp::geqrf( A->rows(), A->cols(),A->col_ptr(1), A->max_rows()
    , tau->raw_ptr(), work->raw_ptr(), work->dim(), &info  );
  if( info != 0 ) {
    TEST_FOR_EXCEPTION(
      info < 0, std::invalid_argument
      ,"geqrf(...): Error, Invalid argument "
      << -info << " sent to LAPACK function xGEQRF(...)" );
  }
  // If you get here A and tau contain the QR factors of input A.
}

void DenseLinAlgLAPack::ormrq(
  BLAS_Cpp::Side side, BLAS_Cpp::Transp trans
  ,const DMatrixSlice& A, const DVectorSlice& tau
  ,DMatrixSlice* C, DVectorSlice* work
  )
{
  FortranTypes::f_int info;
  LAPACK_Cpp::ormqr( side, trans, C->rows(), C->cols()
    , tau.dim(), A.col_ptr(1), A.max_rows()
    , tau.raw_ptr(), C->col_ptr(1), C->max_rows()
    , work->raw_ptr(), work->dim(), &info );
  if( info != 0 ) {
    TEST_FOR_EXCEPTION(
      info < 0, std::invalid_argument
      ,"ormrq(...): Error, Invalid argument "
      << -info << " sent to LAPACK function xORMRQ(...)" );
  }
  // If you get here C contains the desired matrix-matrix multiplication.
}

void DenseLinAlgLAPack::sytrf(
  DMatrixSliceTriEle* A, FortranTypes::f_int ipiv[]
  ,DVectorSlice* work
  )
{
  FortranTypes::f_int info;
  LAPACK_Cpp::sytrf( A->uplo(), A->rows(), A->gms().col_ptr(1)
    , A->gms().max_rows(), ipiv, work->raw_ptr(), work->dim()
    , &info );
  TEST_FOR_EXCEPTION(
    info < 0, std::invalid_argument
    ,"sytrf(...): Error, Invalid argument "
    << -info << " sent to LAPACK function xSYTRF(...)" );
  TEST_FOR_EXCEPTION(
    info > 0, FactorizationException
    ,"sytrf(...): Error, xSYTRF(...) indicates a singular matrix, "
    << "D("<<info<<","<<info<<") is zero." );
  // If we get here A and ipiv contain the factorization.
}

void DenseLinAlgLAPack::sytrs(
  const DMatrixSliceTriEle& A, FortranTypes::f_int ipiv[]
  ,DMatrixSlice* B, DVectorSlice* work
  )
{
  TEST_FOR_EXCEPTION(
    (A.rows() != B->rows()), std::invalid_argument
    ,"sytrs(...) : Error, The number of rows in A and B must match."
    );
  FortranTypes::f_int info;
  LAPACK_Cpp::sytrs(  A.uplo(), A.rows(), B->cols(), A.gms().col_ptr(1)
    , A.gms().max_rows(), ipiv, B->col_ptr(1), B->max_rows()
    , &info );
  TEST_FOR_EXCEPTION(
    info < 0, std::invalid_argument
    ,"sytrs(...): Error, Invalid argument "
    << -info << " sent to LAPACK function xSYTRS(...)"
    );
}

void DenseLinAlgLAPack::getrf(
  DMatrixSlice* A, FortranTypes::f_int ipiv[], FortranTypes::f_int* rank
  )
{
  FortranTypes::f_int info;
  LAPACK_Cpp::getrf(
    A->rows(), A->cols(), A->col_ptr(1), A->max_rows()
    ,ipiv, &info
    );
  *rank = my_min( A->rows(), A->cols() );
  TEST_FOR_EXCEPTION(
    info < 0, std::invalid_argument
    ,"getrf(...): Error, Invalid argument "
    << -info << " sent to LAPACK function xGETRF(...)" );
  if(info > 0)
    *rank = info - 1;
}

void DenseLinAlgLAPack::getrs(
  const DMatrixSlice& LU, const FortranTypes::f_int ipiv[], BLAS_Cpp::Transp transp
  ,DMatrixSlice* B
  )
{
  TEST_FOR_EXCEPTION(
    (LU.rows() != LU.cols() || LU.rows() != B->rows() ), std::invalid_argument
    ,"getrs(...) : Error, A must be square and the number of rows in A and B must match."
    );
  FortranTypes::f_int info;
  LAPACK_Cpp::getrs(
    transp, LU.rows(), B->cols(), LU.col_ptr(1), LU.max_rows(), ipiv
    ,B->col_ptr(1), B->max_rows(), &info
    );
  TEST_FOR_EXCEPTION(
    info < 0, std::invalid_argument
    ,"getrs(...): Error, Invalid argument "
    << -info << " sent to LAPACK function xGETRS(...)" );
  // If we get here B is the solution to the linear system.
}
