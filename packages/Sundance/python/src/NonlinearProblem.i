// -*- c++ -*-


%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceNonlinearProblem.hpp"

  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"


 // SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"

namespace Sundance
{
class NonlinearProblem
{
public:
  /** Empty ctor */
  NonlinearProblem();

  /** Construct with a mesh, equation set, bcs, test and unknown funcs,
   * and a vector type */
  NonlinearProblem(const Sundance::Mesh& mesh, const Sundance::Expr& eqn, const Sundance::Expr& bc,
    const Sundance::Expr& test, const Sundance::Expr& unk, const Sundance::Expr& u0, 
    const TSFExtended::VectorType<double>& vecType);

  /** Construct with a mesh, equation set, bcs, test and unknown funcs,
   * parameters, and a vector type */
  NonlinearProblem(const Sundance::Mesh& mesh, const Sundance::Expr& eqn, const Sundance::Expr& bc,
    const Sundance::Expr& test, const Sundance::Expr& unk, const Sundance::Expr& u0, 
    const Sundance::Expr& params, const Sundance::Expr& paramVals,  
    const TSFExtended::VectorType<double>& vecType);


  

  /** Compute direct sensitivities to parameters */
  Sundance::Expr computeSensitivities(const TSFExtended::LinearSolver<double>& solver) const ;

  /** Solve the nonlinear problem */
  NOX::StatusTest::StatusType solve(const TSFExtended::NOXSolver& solver) const ;

  /** Return the current evaluation point as a Sundance expression */
  Sundance::Expr getU0() const ;

  /** Set an initial guess */
  void setInitialGuess(const Sundance::Expr& u0New);
};

}

