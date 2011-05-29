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

#include "TSFDenseSerialMatrix.hpp"
#include "TSFDenseSerialMatrixFactory.hpp"
#include "TSFSerialVector.hpp"
#include "TSFVectorSpaceDecl.hpp"  
#include "TSFVectorDecl.hpp"
#include "TSFLinearOperatorDecl.hpp"
#include "Teuchos_BLAS.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#include "TSFVectorImpl.hpp"
#endif

extern "C"
{
void dgesv_(int *n, int *nrhs, double *a, int* lda, 
  int *ipiv, double *b, int *ldb, int *info);

void dgesvd_( char* jobu, char* jobvt, int* m, int* n, double* a,
  int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt,
  double* work, int* lwork, int* info );
}
using std::max;
using std::min;

using namespace TSFExtended;
using namespace Teuchos;
using namespace Thyra;
using std::setw;

DenseSerialMatrix::DenseSerialMatrix(
  const RCP<const SerialVectorSpace>& domain,
  const RCP<const SerialVectorSpace>& range)
  : range_(range),
    domain_(domain),
    nRows_(range_->dim()),
    nCols_(domain_->dim()),
    data_(nRows_*nCols_)
{}


Teuchos::ETransp thyraTransToTeuchosTrans(const Thyra::EOpTransp M_trans)
{
  switch(M_trans)
  {
    case Thyra::NOTRANS:
      return Teuchos::NO_TRANS;
    case Thyra::TRANS:
      return Teuchos::TRANS;
    case Thyra::CONJTRANS:
      return Teuchos::CONJ_TRANS;
    default:
      TEST_FOR_EXCEPT(true);
  }
  return Teuchos::NO_TRANS; // -Wall
}

void DenseSerialMatrix::applyOp(
  const Thyra::EOpTransp M_trans,
  const Vector<double>& in,
  Vector<double> out) const
{
  const SerialVector* rvIn = SerialVector::getConcrete(in);
  SerialVector* rvOut = SerialVector::getConcrete(out);

  Teuchos::BLAS<OrdType, double> blas;
  int lda = numRows();
  Teuchos::ETransp trans = thyraTransToTeuchosTrans(M_trans);
  blas.GEMV(trans, numRows(), numCols(), 1.0, dataPtr(), 
    lda, rvIn->dataPtr(), 1, 1.0, rvOut->dataPtr(), 1);
}

void DenseSerialMatrix::addToRow(int globalRowIndex,
  int nElemsToInsert,
  const int* globalColumnIndices,
  const double* elementValues)
{
  int r = globalRowIndex;
  for (int k=0; k<nElemsToInsert; k++)
  {
    int c = globalColumnIndices[k];
    double x = elementValues[k];
    data_[r + c*numRows()] = x;
  }
}

void DenseSerialMatrix::zero()
{
  for (int i=0; i<data_.size(); i++) data_[i] = 0.0;
}


void DenseSerialMatrix::print(std::ostream& os) const
{
  if (numCols() <= 4)
  {
    for (int i=0; i<numRows(); i++)
    {
      for (int j=0; j<numCols(); j++)
      {
        os << setw(16) << data_[i+numRows()*j];
      }
      os << std::endl;
    }
  }
  else
  {
    for (int i=0; i<numRows(); i++)
    {
      for (int j=0; j<numCols(); j++)
      {
        os << setw(6) << i << setw(6) << j << setw(20) << data_[i+numRows()*j]
           << std::endl;
      }
    }
  }
}

void DenseSerialMatrix::setRow(int row, const Array<double>& rowVals)
{
  TEST_FOR_EXCEPT(rowVals.size() != numCols());
  TEST_FOR_EXCEPT(row < 0);
  TEST_FOR_EXCEPT(row >= numRows());

  for (int i=0; i<rowVals.size(); i++)
  {
    data_[row+numRows()*i] = rowVals[i];
  }
}


namespace TSFExtended
{


SolverState<double> denseSolve(const LinearOperator<double>& A,
  const Vector<double>& b,
  Vector<double>& x)
{
  const DenseSerialMatrix* Aptr 
    = dynamic_cast<const DenseSerialMatrix*>(A.ptr().get());
  TEST_FOR_EXCEPT(Aptr==0);
  /* make a working copy, because dgesv will overwrite the matrix */
  DenseSerialMatrix tmp = *Aptr;
  /* Allocate a vector for the solution */
  x = b.copy();
  SerialVector* xptr 
    = dynamic_cast<SerialVector*>(x.ptr().get());
  
  int N = Aptr->numRows();
  int nRHS = 1;
  int LDA = N;
  Array<int> iPiv(N);
  int LDB = N;
  int info = 0;
  dgesv_(&N, &nRHS, tmp.dataPtr(), &LDA, &(iPiv[0]), xptr->dataPtr(),
    &LDB, &info);

  if (info == 0)
  {
    return SolverState<double>(SolveConverged, "solve OK",
      0, 0.0);
  }
  else 
  {
    return SolverState<double>(SolveCrashed, "solve crashed with dgesv info="
      + Teuchos::toString(info),
      0, 0.0);
  }
}


void denseSVD(const LinearOperator<double>& A,
  LinearOperator<double>& U,  
  Vector<double>& Sigma,
  LinearOperator<double>& Vt)
{
  VectorSpace<double> mSpace = A.range();
  RCP<const SerialVectorSpace> rmSpace 
    = rcp_dynamic_cast<const SerialVectorSpace>(mSpace.ptr());

  VectorSpace<double> nSpace = A.domain();
  RCP<const SerialVectorSpace> rnSpace 
    = rcp_dynamic_cast<const SerialVectorSpace>(nSpace.ptr());

  const DenseSerialMatrix* Aptr 
    = dynamic_cast<const DenseSerialMatrix*>(A.ptr().get());
  TEST_FOR_EXCEPT(Aptr==0);
  /* make a working copy, because dgesvd will overwrite the matrix */
  DenseSerialMatrix ATmp = *Aptr;

  int M = ATmp.numRows();
  int N = ATmp.numCols();
  int S = min(M, N);
  
  VectorSpace<double> sSpace;
  if (S==M) sSpace = mSpace;
  else sSpace = nSpace;

  RCP<const SerialVectorSpace> rsSpace 
    = rcp_dynamic_cast<const SerialVectorSpace>(sSpace.ptr());

  Sigma = sSpace.createMember();
  SerialVector* sigPtr
    = dynamic_cast<SerialVector*>(Sigma.ptr().get());
  TEST_FOR_EXCEPT(sigPtr==0);

  DenseSerialMatrixFactory umf(rsSpace, rmSpace);
  DenseSerialMatrixFactory vtmf(rnSpace, rsSpace);
  
  U = umf.createMatrix();
  Vt = vtmf.createMatrix();

  DenseSerialMatrix* UPtr 
    = dynamic_cast<DenseSerialMatrix*>(U.ptr().get());
  TEST_FOR_EXCEPT(UPtr==0);

  DenseSerialMatrix* VtPtr 
    = dynamic_cast<DenseSerialMatrix*>(Vt.ptr().get());
  TEST_FOR_EXCEPT(VtPtr==0);
  
  double* uData = UPtr->dataPtr();
  double* vtData = VtPtr->dataPtr();
  double* aData = ATmp.dataPtr();
  double* sData = sigPtr->dataPtr();

  char jobu = 'S';
  char jobvt = 'S';
 
  int LDA = M;
  int LDU = M;
  int LDVT = S;

  int LWORK = max(1, max(3*min(M,N)+max(M,N), 5*min(M,N)));
  Array<double> work(LWORK);
  
  int info = 0;

  dgesvd_(&jobu, &jobvt, &M, &N, aData, &LDA, sData, uData, &LDU, 
    vtData, &LDVT, &(work[0]), &LWORK, &info);

  TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
    "dgesvd failed with error code info=" << info);

  
  
}

}
