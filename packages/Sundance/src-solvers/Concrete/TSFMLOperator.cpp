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

#include "TSFMLOperator.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_MPIComm.hpp"
#include "TSFEpetraVector.hpp"



#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorImpl.hpp"
#include "TSFLinearOperatorImpl.hpp"
#endif

using namespace TSFExtended;
using namespace Teuchos;

MLOperator::MLOperator(
  const LinearOperator<double>& op,
  const ParameterList& mlParams)
  : mlPrec_(),
    domain_(op.domain().ptr()),
    range_(op.range().ptr())
{
	Epetra_CrsMatrix& A = EpetraMatrix::getConcrete(op);
  
  
  mlPrec_ = rcp(new ML_Epetra::MultiLevelPreconditioner(A, mlParams));
}


void MLOperator::generalApply(const Thyra::EOpTransp M_trans,
  const Thyra::VectorBase<double>& x,
  Thyra::VectorBase<double>* y,
  const double alpha,
  const double beta) const
{
  /* grab the epetra vector objects underlying the input and output vectors */
  const EpetraVector* epIn = dynamic_cast<const EpetraVector*>(&x);
  TEST_FOR_EXCEPTION(epIn == 0, std::runtime_error,
                     "IfpackOperator apply: input vector is "
                     "not an EpetraVector");

  const Epetra_Vector* in = epIn->epetraVec().get();

  EpetraVector* epy = dynamic_cast<EpetraVector*>(y);
  TEST_FOR_EXCEPTION(epy == 0, std::runtime_error,
                     "IfpackOperator apply: output vector is "
                     "not an EpetraVector");

  Epetra_Vector* yy = epy->epetraVec().get();

  /* if beta != 0, we have to create a temporary vector to hold the
   * intermediate result beta*y. */
  RCP<Thyra::VectorBase<double> > tsfOut;
  Epetra_Vector* tmp;

  if (beta!=0.0) /* we need to have storage for the result of op*x */
    {
      tsfOut = createMember(y->space());
      EpetraVector* ep = dynamic_cast<EpetraVector*>(tsfOut.get());
      tmp = ep->epetraVec().get();
    }
  else /* we overwrite y with the application of the op */
    {
      tmp = yy;
    }


  int ierr;

  /* do the solve (or transpose solve) */
  if (M_trans==NOTRANS)
    {
      ierr = mlPrec_->ApplyInverse(*in, *tmp);
    }
  else
    {
      TEST_FOR_EXCEPTION(M_trans != NOTRANS, std::runtime_error,
        "ML preconditioner does not support transposes");
    }

  /* if necessary, add beta*y */
  if (beta != 0.0)
    {
      /* Compute yy = alpha*tmp + beta*yy */
      yy->Update(alpha, *tmp, beta);
    }  
  else if (alpha != 1.0) /* compute yy = alpha*tmp */
    {
      /* Because beta != 0.0, an earlier conditional has set tmp=yy. 
       * We can therefore do the computation on tmp and it will 
       * modify the contents of yy as a side effect. */
      tmp->Scale(alpha);
    }

  /* At this point, the contents of y should be yy. We are done. */
}
  
