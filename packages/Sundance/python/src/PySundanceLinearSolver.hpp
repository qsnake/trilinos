#ifndef PYSUNDANCELINEARSOLVER_H
#define PYSUNDANCELINEARSOLVER_H

#include "PySundanceCallback.hpp"
#include "TSFLinearSolverBaseDecl.hpp"
#include "SundanceHandleable.hpp"


namespace TSFExtended
{
  class PySundanceLinearSolver : public LinearSolverBase<double>,
                                 public Sundance::Handleable<LinearSolverBase<double> >
  {
  public:
    /** */
    PySundanceLinearSolver(PyObject* functor);
    /** */
    ~PySundanceLinearSolver();

    /* */
    GET_RCP(LinearSolverBase<double>);

    /** */
    SolverState<double> solve(const LinearOperator<double>& op,
                              const Vector<double>& rhs,
                              Vector<double>& soln) const ;

    std::string description() const {return "I'm a PySundanceLinearSolver";}

    PyObject* pySolve(PyObject* op, PyObject* rhs, PyObject* x0) const ;
  protected:
    PyObject * setSolve(PyObject *);



  private:
    // Private and not implemented
    //PyInterface();
    PySundanceLinearSolver(const PySundanceLinearSolver &);
    PySundanceLinearSolver & operator=(const PySundanceLinearSolver &);

  private:
    PyObject* py_functor_;
    mutable PySundanceCallback  py_solve_;
  };
}

#endif // 
