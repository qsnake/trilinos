#include "PySundanceLinearSolver.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#include "TSFLinearSolverImpl.hpp"
#endif


using namespace TSFExtended;


namespace TSFExtended
{
  SolverState<double> 
  PySundanceLinearSolver_solve(const PySundanceLinearSolver* solver,
                               const LinearOperator<double>& op,
                               const Vector<double>& rhs,
                               Vector<double>& soln);
}

PySundanceLinearSolver::PySundanceLinearSolver(PyObject* functor) 
  : LinearSolverBase<double>(ParameterList()), 
    py_functor_(functor), py_solve_()

{
  // Increment the reference count
  Py_XINCREF(py_functor_);

  // If the python object has a "solve" attribute, set it
  // to be the PySundanceLinearSolver solve() callback function
  if (PyObject_HasAttrString (py_functor_,
			      "solve")) {
    setSolve(PyObject_GetAttrString(py_functor_,
                                    "solve"));
  }
}

PySundanceLinearSolver::~PySundanceLinearSolver() {
  // Decrement the reference count
  Py_XDECREF(py_functor_);
}


SolverState<double> PySundanceLinearSolver::solve(const LinearOperator<double>& op,
                                                  const Vector<double>& rhs,
                                                  Vector<double>& soln) const
{
  return PySundanceLinearSolver_solve(this, op, rhs, soln);
}


PyObject * PySundanceLinearSolver::setSolve(PyObject * p_pyObject)
{
  return py_solve_.setFunction(p_pyObject);
}

PyObject* PySundanceLinearSolver::pySolve(PyObject* opObj, PyObject* rhsObj, 
                                          PyObject* solnObj) const 
{
  PyObject* arglist = Py_BuildValue("(OOO)", opObj, rhsObj, solnObj);

  PyObject* result = PyEval_CallObject(py_solve_.getFunction(), arglist);
  
  Py_DECREF(arglist);  // All done with argument list

  return result;
}
