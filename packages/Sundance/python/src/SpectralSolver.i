// -*- c++ -*-


%{

#define HAVE_PY_FIAT
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceStochBlockJacobiSolver.hpp"

  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"

namespace Sundance
{
  class StochBlockJacobiSolver
  {
  public:
    StochBlockJacobiSolver(
    const TSFExtended::LinearSolver<double>& diagonalSolver,
    const SpectralBasis& pcBasis, 
    double convTol,
    int maxIters,
    int verbosity);

    void solve(const Teuchos::Array<TSFExtended::LinearOperator<double> >& KBlock,
      const Teuchos::Array<int>& hasNonzeroMatrixBlock,
      const Teuchos::Array<TSFExtended::Vector<double> >& fBlock,
      Teuchos::Array<TSFExtended::Vector<double> >& xBlock) const ;

    
    void solve(const Teuchos::Array<TSFExtended::LinearOperator<double> >& KBlock,
      const Teuchos::Array<TSFExtended::Vector<double> >& fBlock,
      Teuchos::Array<TSFExtended::Vector<double> >& xBlock) const ;
  };

}





%template(IntArray) Teuchos::Array<int>;
%template(LinearOperatorArray) Teuchos::Array<TSFExtended::LinearOperator<double> >;



