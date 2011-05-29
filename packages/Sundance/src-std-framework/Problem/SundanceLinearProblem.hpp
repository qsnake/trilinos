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

#ifndef SUNDANCE_LINEARPROBLEM_H
#define SUNDANCE_LINEARPROBLEM_H

#include "SundanceDefs.hpp"
#include "SundanceLinearSolveDriver.hpp"
#include "SundanceObjectWithVerbosity.hpp"

namespace Sundance
{
using namespace Teuchos;

class Assembler;

/** 
 * LinearProblem encapsulates all information needed to form
 * a discrete linear problem. 
 */
class LinearProblem 
{
public:
  /** Empty ctor */
  LinearProblem();
    
  /** Construct with a mesh, equation set, bcs, test and unknown funcs,
   * and a vector type. */
  LinearProblem(const Mesh& mesh, const Expr& eqn, const Expr& bc,
    const Expr& test, const Expr& unk, 
    const TSFExtended::VectorType<double>& vecType
    );
    
  /** Construct with a mesh, equation set, bcs, and blocks of variables */
  LinearProblem(const Mesh& mesh, const Expr& eqn, const Expr& bc,
    const BlockArray& test, const BlockArray& unk);
    
  /** Construct with a mesh, equation set, bcs, test and unknown funcs,
   * parameters, and a vector type. */
  LinearProblem(const Mesh& mesh, const Expr& eqn, const Expr& bc,
    const Expr& test, const Expr& unk, 
    const Expr& unkParams, const Expr& unkParamVals, 
    const TSFExtended::VectorType<double>& vecType);
    
  /** Construct with a mesh, equation set, bcs, parameters, and blocks of
      variables */
  LinearProblem(const Mesh& mesh, const Expr& eqn, const Expr& bc,
    const BlockArray& test, const BlockArray& unk, 
    const Expr& unkParams, const Expr& unkParamVals);

  /** */
  LinearProblem(const RCP<Assembler>& assembler);

  /** Solve the problem using the specified solver algorithm */
  Expr solve(const LinearSolver<double>& solver) const ;

  /** Solve the problem, writing the solution into the given function */
  SolverState<double> solve(const LinearSolver<double>& solver,
    Expr& soln) const ;

  /** Return the multivector on the right-hand side of the linear equation */
  Array<Vector<double> > getRHS() const ;

  /** Return the vector on the right-hand side of the linear equation */
  Vector<double> getSingleRHS() const {return getRHS()[0];}

  /** Return the operator on the left-hand side of the equation */
  LinearOperator<double> getOperator() const ;

  /** Return the map from cells and functions to row indices */
  const RCP<DOFMapBase>& rowMap(int blockRow) const ;
    
  /** Return the map from cells and functions to column indices */
  const RCP<DOFMapBase>& colMap(int blockCol) const ;

  /** Return the discrete space in which solutions live */
  const Array<RCP<DiscreteSpace> >& solnSpace() const ;

    
  /** Return the set of row indices marked as 
   * essential boundary conditions */
  const RCP<Set<int> >& bcRows(int blockRow) const ;

  /** Return the number of block rows in the problem  */
  int numBlockRows() const ;

  /** Return the number of block cols in the problem  */
  int numBlockCols() const ;


  /** Convert from a BC-partitioned solution vector to a 
   * monolithic vector */
  Vector<double> 
  convertToMonolithicVector(const Array<Vector<double> >& internalBlock,
    const Array<Vector<double> >& bcBlock) const ;

  /** */
  Expr formSolutionExpr(const Array<Vector<double> >& vec) const ;

  /** Flag indicating whether to stop on a solve failure */
  static bool& solveFailureIsFatal()
    {return LinearSolveDriver::solveFailureIsFatal();}
    

  /** Flag indicating whether to write out the matrix and vector
   * after a solve failure */
  static bool& dumpBadMatrix() 
    {return LinearSolveDriver::dumpBadMatrix();}

  /** Filename for dump of bad matrix */
  static std::string& badMatrixFilename() 
    {return LinearSolveDriver::badMatrixFilename();}

  /** Filename for dump of bad vector */
  static std::string& badVectorFilename() 
    {return LinearSolveDriver::badVectorFilename();}

    

private:

      
  /** */
  RCP<Assembler> assembler_;

  /** */
  mutable LinearOperator<double> A_;

  /** */
  mutable Array<Vector<double> > rhs_;

  /** */
  Array<Array<string> > names_;

  /** */
  LinearSolveDriver solveDriver_;
    
};
}


#endif
