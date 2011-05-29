#ifndef SUNDANCE_STOCHBLOCKJACOBISOLVER_HPP
#define SUNDANCE_STOCHBLOCKJACOBISOLVER_HPP

#include "TSFLinearSolverDecl.hpp"
#include "SundanceSpectralBasis.hpp"

using TSFExtended::LinearSolver;
using TSFExtended::LinearOperator;
using TSFExtended::Vector;
using Sundance::SpectralBasis;

namespace Sundance
{

class StochBlockJacobiSolver
{
public:
  /** */
  StochBlockJacobiSolver(
    const LinearSolver<double>& diagonalSolver,
    const SpectralBasis& pcBasis, 
    double convTol,
    int maxIters,
    int verbosity)
    : diagonalSolver_(diagonalSolver),
      pcBasis_(pcBasis),
      convTol_(convTol),
      maxIters_(maxIters),
      verbosity_(verbosity)
    {}

  /** */
  void solve(const Array<LinearOperator<double> >& KBlock,
    const Array<int>& hasNonzeroMatrixBlock,
    const Array<Vector<double> >& fBlock,
    Array<Vector<double> >& xBlock) const ;

  /** */
  void solve(const Array<LinearOperator<double> >& KBlock,
    const Array<Vector<double> >& fBlock,
    Array<Vector<double> >& xBlock) const ;

private:
  LinearSolver<double> diagonalSolver_;
  SpectralBasis pcBasis_;
  double convTol_;
  int maxIters_;
  int verbosity_;
};
}

#endif
