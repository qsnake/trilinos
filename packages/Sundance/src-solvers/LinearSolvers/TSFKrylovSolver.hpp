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

#ifndef TSFKRYLOVSOLVER_HPP
#define TSFKRYLOVSOLVER_HPP

#include "SundanceDefs.hpp"
#include "TSFIterativeSolver.hpp"
#include "TSFPreconditionerFactory.hpp"
#include "TSFILUKPreconditionerFactory.hpp"

namespace TSFExtended
{
using namespace Teuchos;

/**
 *
 */
template <class Scalar>
class KrylovSolver : public IterativeSolver<Scalar>
{
public:
  /** */
  KrylovSolver(const ParameterList& params);
  /** */
  KrylovSolver(const ParameterList& params,
    const PreconditionerFactory<Scalar>& precond);

  /** */
  virtual ~KrylovSolver(){;}

  /** */
  virtual SolverState<Scalar> solve(const LinearOperator<Scalar>& op,
    const Vector<Scalar>& rhs,
    Vector<Scalar>& soln) const ;
protected:
  virtual SolverState<Scalar> solveUnprec(const LinearOperator<Scalar>& op,
    const Vector<Scalar>& rhs,
    Vector<Scalar>& soln) const = 0 ;

  const PreconditionerFactory<Scalar>& precond() const {return precond_;}

private:
  PreconditionerFactory<Scalar> precond_;
};

  
template <class Scalar> inline
KrylovSolver<Scalar>::KrylovSolver(const ParameterList& params)
  : IterativeSolver<Scalar>(params), precond_()
{
  if (!params.isParameter("Precond")) return;

  const std::string& precondType = params.template get<string>("Precond");

  if (precondType=="ILUK")
  {
    precond_ = new ILUKPreconditionerFactory<Scalar>(params);
  }
}

template <class Scalar> inline
KrylovSolver<Scalar>::KrylovSolver(const ParameterList& params,
  const PreconditionerFactory<Scalar>& precond)
  : IterativeSolver<Scalar>(params), precond_(precond)
{
  TEST_FOR_EXCEPTION(params.isParameter("Precond"), std::runtime_error,
    "ambiguous preconditioner specification in "
    "KrylovSolver ctor: parameters specify "
    << params.template get<string>("Precond") 
    << " but preconditioner argument is " 
    << precond);
}

template <class Scalar> inline
SolverState<Scalar> KrylovSolver<Scalar>
::solve(const LinearOperator<Scalar>& op,
  const Vector<Scalar>& rhs,
  Vector<Scalar>& soln) const
{
  if (precond_.ptr().get()==0) 
  {
    return solveUnprec(op, rhs, soln);
  }


  Preconditioner<Scalar> p = precond_.createPreconditioner(op);
    
  if (!p.hasRight())
  {
    LinearOperator<Scalar> A = p.left()*op;
    Vector<Scalar> newRHS = rhs.space().createMember();
    p.left().apply(rhs, newRHS);
    return solveUnprec(A, newRHS, soln);
  }
  else if (!p.hasLeft())
  {
    LinearOperator<Scalar> A = op * p.right();
    Vector<Scalar> intermediateSoln;
    SolverState<Scalar> rtn 
      = solveUnprec(A, rhs, intermediateSoln);
    if (rtn.finalState()==SolveConverged) 
    {
      p.right().apply(intermediateSoln, soln);
    }
    return rtn;
  }
  else
  {
    LinearOperator<Scalar> A = p.left() * op * p.right();
    Vector<Scalar> newRHS;
    p.left().apply(rhs, newRHS);
    Vector<Scalar> intermediateSoln;
    SolverState<Scalar> rtn 
      = solveUnprec(A, newRHS, intermediateSoln);
    if (rtn.finalState()==SolveConverged) 
    {
      p.right().apply(intermediateSoln, soln);
    }
    return rtn;
  }
}
  
}

#endif

