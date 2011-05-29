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
#include "Teuchos_Array.hpp"
#include "Teuchos_MPIComm.hpp"
#include "TSFIfpackOperator.hpp"
#include "TSFGenericLeftPreconditioner.hpp"
#include "TSFGenericRightPreconditioner.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_getConst.hpp"
#include "EpetraTSFOperator.hpp"
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"
#include "Thyra_Config.h"

#ifdef HAVE_THYRA_EPETRA
#include "Thyra_EpetraThyraWrappers.hpp"
#endif

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorImpl.hpp"
#include "TSFLinearCombinationImpl.hpp"
#include "TSFLinearOperatorImpl.hpp"
#include "TSFLinearSolverImpl.hpp"
#endif

using namespace TSFExtended;
using namespace Teuchos;
using namespace Thyra;
using namespace Epetra;

Epetra_TSFOperator::Epetra_TSFOperator(const LinearOperator<double>& A,
				       const LinearSolver<double>& solver)
  : A_(A), solver_(solver), useTranspose_(false), comm_(), domain_(), range_(),
    isNativeEpetra_(false), isCompoundEpetra_(false), label_(A.description())
{
  const EpetraMatrix* em = dynamic_cast<const EpetraMatrix*>(A.ptr().get());
  const EpetraVectorSpace* ed = dynamic_cast<const EpetraVectorSpace*>(A.domain().ptr().get());
  const EpetraVectorSpace* er = dynamic_cast<const EpetraVectorSpace*>(A.range().ptr().get());

  if (em)
    {
      isNativeEpetra_ = true;
      const Epetra_CrsMatrix* crs = em->crsMatrix();
      domain_ = rcp(new Epetra_Map(crs->OperatorDomainMap()));
      range_ = rcp(new Epetra_Map(crs->OperatorRangeMap()));
      useTranspose_ = crs->UseTranspose();
      comm_ = rcp(crs->Comm().Clone());
    }
  else if (er != 0 && ed != 0)
    {
      domain_ = ed->epetraMap();
      range_ = er->epetraMap();
      comm_ = rcp(domain_->Comm().Clone());
      isCompoundEpetra_ = true;
    }
  else
    {
      TEST_FOR_EXCEPT(true);
    }
}


int Epetra_TSFOperator::Apply(const Epetra_MultiVector& in, Epetra_MultiVector& out) const
{
  if (isNativeEpetra_)
    {
      const EpetraMatrix* em = dynamic_cast<const EpetraMatrix*>(A_.ptr().get());
      return em->crsMatrix()->Multiply(useTranspose_, in, out);
    }
  else if (isCompoundEpetra_)
    {
      const Epetra_Vector* cevIn = dynamic_cast<const Epetra_Vector*>(&in);
      Epetra_Vector* evIn = const_cast<Epetra_Vector*>(cevIn);
      Epetra_Vector* evOut = dynamic_cast<Epetra_Vector*>(&out);
      TEST_FOR_EXCEPTION(evIn==0, std::runtime_error, "Epetra_TSFOperator::Apply "
			 "cannot deal with multivectors");
      TEST_FOR_EXCEPTION(evOut==0, std::runtime_error, "Epetra_TSFOperator::Apply "
			 "cannot deal with multivectors");

      const EpetraVectorSpace* ed 
	= dynamic_cast<const EpetraVectorSpace*>(A_.domain().ptr().get());
      const EpetraVectorSpace* er 
	= dynamic_cast<const EpetraVectorSpace*>(A_.range().ptr().get());
      TEST_FOR_EXCEPTION(er == 0 || ed==0, std::runtime_error, 
			 "this should never happen, because we have found "
			 "Epetra domain and range in the ctor");

      RCP<Thyra::VectorBase<double> > vpIn 
	= rcp(new EpetraVector(rcp(ed, false), 
			       rcp(evIn, false)));
      RCP<Thyra::VectorBase<double> > vpOut 
	= rcp(new EpetraVector(rcp(er, false), 
			       rcp(evOut, false)));
      Vector<double> vIn = vpIn;
      Vector<double> vOut = vpOut;

      A_.apply(vIn, vOut);
      out = EpetraVector::getConcrete(vOut);
      return 0;
    }
  else
    {
      TEST_FOR_EXCEPT(true);
      return -1; // -Wall
    }
}

int Epetra_TSFOperator::ApplyInverse(const Epetra_MultiVector& in, Epetra_MultiVector& out) const
{
  
  TEST_FOR_EXCEPTION(solver_.ptr().get()==0, std::runtime_error,
		     "no solver provided for Epetra_TSFOperator::ApplyInverse");
  TEST_FOR_EXCEPTION(!isNativeEpetra_ && !isCompoundEpetra_, std::runtime_error,
		     "Epetra_TSFOperator::ApplyInverse expects either "
		     "a native epetra operator or a compound operator with "
		     "Epetra domain and range spaces");
  const Epetra_Vector* cevIn = dynamic_cast<const Epetra_Vector*>(&in);
  Epetra_Vector* evIn = const_cast<Epetra_Vector*>(cevIn);
  Epetra_Vector* evOut = dynamic_cast<Epetra_Vector*>(&out);

  TEST_FOR_EXCEPTION(evIn==0, std::runtime_error, "Epetra_TSFOperator::Apply "
		     "cannot deal with multivectors");
  TEST_FOR_EXCEPTION(evOut==0, std::runtime_error, "Epetra_TSFOperator::Apply "
		     "cannot deal with multivectors");

  const EpetraVectorSpace* ed 
    = dynamic_cast<const EpetraVectorSpace*>(A_.range().ptr().get());
  const EpetraVectorSpace* er 
    = dynamic_cast<const EpetraVectorSpace*>(A_.domain().ptr().get());

  RCP<Thyra::VectorBase<double> > vpIn 
    = rcp(new EpetraVector(rcp(ed, false), 
			   rcp(evIn, false)));
  RCP<Thyra::VectorBase<double> > vpOut 
    = rcp(new EpetraVector(rcp(er, false), 
			   rcp(evOut, false)));
  Vector<double> vIn = vpIn;
  Vector<double> vOut = vpOut;
  
  SolverState<double> state = solver_.solve(A_, vIn, vOut);

  if (state.finalState() == SolveCrashed) return -1;
  else if (state.finalState() == SolveFailedToConverge) return -2;
  else out = EpetraVector::getConcrete(vOut);

  return 0;
}




double Epetra_TSFOperator::NormInf() const 
{
  TEST_FOR_EXCEPT(true);
  return -1; // -Wall
}

const char* Epetra_TSFOperator::Label() const 
{
  return label_.c_str();
}
