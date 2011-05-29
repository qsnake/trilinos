/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
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
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#ifndef SUNDANCE_NONLINEARPROBLEM_H
#define SUNDANCE_NONLINEARPROBLEM_H

#include "SundanceDefs.hpp"
#include "SundanceNLOp.hpp"
#include "TSFNOXSolver.H"

namespace Sundance
{
using namespace Teuchos;


/** 
 * NonlinearProblem encapsulates a discrete nonlinear problem
 */
class NonlinearProblem 
  : public ObjectWithClassVerbosity<NonlinearProblem>
{
public:
  /** Empty ctor */
  NonlinearProblem();

  /** Construct with a mesh, equation set, bcs, test and unknown funcs,
   * and a vector type */
  NonlinearProblem(const Mesh& mesh, const Expr& eqn, const Expr& bc,
    const Expr& test, const Expr& unk, const Expr& u0, 
    const TSFExtended::VectorType<double>& vecType);

  /** Construct with a mesh, equation set, bcs, test and unknown funcs,
   * parameters, and a vector type */
  NonlinearProblem(const Mesh& mesh, const Expr& eqn, const Expr& bc,
    const Expr& test, const Expr& unk, const Expr& u0, 
    const Expr& params, const Expr& paramVals,  
    const TSFExtended::VectorType<double>& vecType);


  /** */
  NonlinearProblem(const RCP<Assembler>& assembler, 
    const Expr& u0);

  /** Compute direct sensitivities to parameters */
  Expr computeSensitivities(const LinearSolver<double>& solver) const 
    {return op_->computeSensitivities(solver);}

  /** Solve the nonlinear problem */
  NOX::StatusTest::StatusType solve(const NOXSolver& solver) const ;

  /** Return the current evaluation point as a Sundance expression */
  Expr getU0() const {return op_->getU0();}

  /** Set an initial guess */
  void setInitialGuess(const Expr& u0New) {op_->setInitialGuess(u0New);}
      

  /** Compute the residual and Jacobian at the current evaluation point */
  LinearOperator<double> computeJacobianAndFunction(Vector<double>& functionValue) const 
    {return op_->computeJacobianAndFunction(functionValue);}
      
  /** Write the Jacobian and residual into the objects provided */
  void computeJacobianAndFunction(LinearOperator<double>& J,
    Vector<double>& resid) const 
    {op_->computeJacobianAndFunction(J, resid);}

  /** Compute the residual at the current eval point */
  TSFExtended::Vector<double> computeFunctionValue() const 
    {return op_->computeFunctionValue();}
      
  /** Write the residual into the object provided */
  void computeFunctionValue(Vector<double>& resid) const 
    {op_->computeFunctionValue(resid);}
      
  /** Get an initial guess */
  TSFExtended::Vector<double> getInitialGuess() const 
    {return op_->getInitialGuess();}
      
  /** Create the Jacobian object, but don't fill it in. */
  LinearOperator<double> allocateJacobian() const 
    {return op_->allocateJacobian();}

private:
  RCP<NLOp> op_;
};
}




#endif
