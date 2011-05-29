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

#ifndef TSFLINEARSOLVERBASEIMPL_HPP
#define TSFLINEARSOLVERBASEIMPL_HPP

#include "TSFLinearSolverBaseDecl.hpp"
#include "TSFPreconditioner.hpp"
#include "TSFPreconditionerFactory.hpp"
#include "Teuchos_ParameterList.hpp"


using namespace Teuchos;


namespace TSFExtended
{

template <class Scalar> inline
const ParameterList& LinearSolverBase<Scalar>::parameters() const 
{return params_;}

template <class Scalar> inline
LinearSolverBase<Scalar>::LinearSolverBase(const ParameterList& params)
  : DefaultObjectWithVerbosity(), 
    params_(params) 
{}

template <class Scalar> inline
ParameterList& LinearSolverBase<Scalar>::parameters() {return params_;}


template <class Scalar> inline
string LinearSolverBase<Scalar>::verbosityParam() {return "Verbosity";}

template <class Scalar>
template <typename T> inline
void LinearSolverBase<Scalar>::setParameter(const ParameterList& params,
  T* dataPtr,
  const std::string& name)
{
  if (!params.isParameter(name)) return;

  TEST_FOR_EXCEPTION(!params.template isType<T>(name), std::runtime_error,
    "invalid type for parameter [" << name << "]"); 

  *dataPtr = params.template get<T>(name);
}

template <class Scalar> inline
void LinearSolverBase<Scalar>::setUserPrec(const PreconditionerFactory<Scalar>& pf)
{
  TEST_FOR_EXCEPTION(true, std::runtime_error,
    "User-defined preconditioning not allowed for generic "
    "linear solver subtypes");
}

template <class Scalar> inline
void LinearSolverBase<Scalar>::setUserPrec(const LinearOperator<Scalar>& P,
  const LinearSolver<Scalar>& pSolver)
{
  TEST_FOR_EXCEPTION(true, std::runtime_error,
    "User-defined preconditioning not allowed for generic "
    "linear solver subtypes");
}

}

#endif
