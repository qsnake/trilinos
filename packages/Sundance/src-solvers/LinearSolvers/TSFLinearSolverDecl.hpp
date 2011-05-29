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

#ifndef TSFLINEARSOLVERDECL_HPP
#define TSFLINEARSOLVERDECL_HPP

#include "SundanceTabs.hpp"
#include "SundanceHandle.hpp"
#include "SundanceHandleable.hpp"
#include "TSFLinearSolverBaseDecl.hpp"
#include "Teuchos_TimeMonitor.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearSolverBaseImpl.hpp"
#include "TSFLinearOperatorImpl.hpp"
#endif

inline static Teuchos::Time& solveTimer() 
{
  static Teuchos::RCP<Teuchos::Time> rtn 
    = Teuchos::TimeMonitor::getNewTimer("linear solve"); 
  return *rtn;
}

namespace TSFExtended
{
using namespace Teuchos;
using namespace Sundance;
  
/**
 *
 */
template <class Scalar>
class LinearSolver : public Sundance::Handle<LinearSolverBase<Scalar> >
{
public:
  /** */
  LinearSolver() : Sundance::Handle<LinearSolverBase<Scalar> >() {;}
  /** */
  LinearSolver( Sundance::Handleable<LinearSolverBase<Scalar> >* rawPtr) 
    : Sundance::Handle<LinearSolverBase<Scalar> >(rawPtr) {;}
  /** */
  LinearSolver(const RCP<LinearSolverBase<Scalar> >& smartPtr)
    : Sundance::Handle<LinearSolverBase<Scalar> >(smartPtr) {;}


  /** Change the convergence tolerance. Default does nothing. */
  void updateTolerance(const double& tol) {this->ptr()->updateTolerance(tol);}

  /** Set a user-defined preconditioner */
  void setUserPrec(const LinearOperator<Scalar>& op,
    const LinearSolver<Scalar>& pSolver) ;

  /** Set a user-defined preconditioner */
  void setUserPrec(const PreconditionerFactory<Scalar>& pf);


  /** */
  SolverState<Scalar> solve(const LinearOperator<Scalar>& op,
    const Vector<Scalar>& rhs,
    Vector<Scalar>& soln) const ;
    
    

  /** */
  const ParameterList& parameters() const ;

  /** */
  ParameterList& parameters() ;
};

  
template <class Scalar> inline 
SolverState<Scalar> LinearSolver<Scalar>
::solve(const LinearOperator<Scalar>& op,
  const Vector<Scalar>& rhs,
  Vector<Scalar>& soln) const
{
  Tabs tab;
  TEST_FOR_EXCEPTION(this->ptr().get()==0, std::runtime_error,
    "null pointer in LinearSolver<Scalar>::solve()");

  TEST_FOR_EXCEPTION(rhs.ptr().get()==0, std::runtime_error,
    "null rhs pointer in LinearSolver<Scalar>::solve()");

  TEST_FOR_EXCEPTION(op.ptr().get()==0, std::runtime_error,
    "null op pointer in LinearSolver<Scalar>::solve()");

  TimeMonitor timer(solveTimer());

  SUNDANCE_MSG1(this->ptr()->verb(), 
    tab << "Solver(" << this->description() << ") starting solve");

  SolverState<Scalar> rtn = this->ptr()->solve(op, rhs, soln);

  SUNDANCE_MSG1(this->ptr()->verb(), 
    tab << "Solver(" << this->description() << ") done solve:");
  Tabs tab1;
  SUNDANCE_MSG2(this->ptr()->verb(), 
    tab << "state=" << rtn);

  return rtn;    
}

template <class Scalar> inline 
const ParameterList& LinearSolver<Scalar>::parameters() const 
{
  TEST_FOR_EXCEPTION(this->ptr().get()==0, std::runtime_error,
    "null pointer in LinearSolver<Scalar>::parameters()");
  return this->ptr()->parameters();
}

template <class Scalar> inline 
ParameterList& LinearSolver<Scalar>::parameters() 
{
  TEST_FOR_EXCEPTION(this->ptr().get()==0, std::runtime_error,
    "null pointer in LinearSolver<Scalar>::parameters()");
  return this->ptr()->parameters();
}

  

  

}

#endif
