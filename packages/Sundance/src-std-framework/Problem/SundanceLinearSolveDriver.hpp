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

#ifndef SUNDANCE_INTERNALSOLVEMANAGER_HPP
#define SUNDANCE_INTERNALSOLVEMANAGER_HPP

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceExpr.hpp"
#include "SundanceBlock.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "TSFLinearOperatorDecl.hpp"
#include "TSFVectorDecl.hpp"
#include "TSFVectorType.hpp"
#include "TSFSolverState.hpp"

namespace Sundance
{

class LinearSolveDriver
{
public:
  /** */
  LinearSolveDriver() {}
    

  /** */
  Expr formSolutionExpr(const Array<Vector<double> >& solnVector,
    const Array<RCP<DiscreteSpace> >& solutionSpace,
    const Array<Array<string> >& names,
    int verb) const ;

  /** */
  SolverState<double> solve(const LinearSolver<double>& solver,
    const LinearOperator<double>& A,
    const Array<Vector<double> >& rhs,
    const Array<RCP<DiscreteSpace> >& solutionSpace,
    const Array<Array<string> >& names,
    int verb,
    Expr& soln) const ;


  /** */
  void writeIntoSolutionExpr(
    const Array<Vector<double> >& solnVector,
    Expr soln, int verb) const ;

  /** Filename for dump of bad matrix */
  static std::string& badMatrixFilename() 
    {static std::string rtn = "badMatrix.dat"; return rtn;}

  /** Filename for dump of bad vector */
  static std::string& badVectorFilename() 
    {static std::string rtn = "badVector.dat"; return rtn;}

  /** Whether a solve failure throws an exception */
  static bool& solveFailureIsFatal()
    {static bool rtn=true; return rtn;}

  /** Whether to dump a matrix upon solve failure */
  static bool& dumpBadMatrix()
    {static bool rtn=false; return rtn;}

  
  

private:
};

}

#endif
