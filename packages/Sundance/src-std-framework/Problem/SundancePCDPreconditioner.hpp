#ifndef SUNDANCE_PCD_PRECONDITIONER_HPP
#define SUNDANCE_PCD_PRECONDITIONER_HPP

#include "SundanceLinearProblem.hpp"
#include "TSFLinearSolverDecl.hpp"
#include "TSFPreconditionerFactory.hpp"


namespace TSFExtended
{
using namespace Sundance;

class PCDPreconditionerFactory
  : public PreconditionerFactoryBase<double>
{
public:
  /** */
  PCDPreconditionerFactory(
    const ParameterList& params,
    const LinearProblem& MpProb,
    const LinearProblem& ApProb,
    const LinearProblem& FpProb
    );


  /** */
  Preconditioner<double> 
  createPreconditioner(const LinearOperator<double>& A) const ;
  /* */
  GET_RCP( PreconditionerFactoryBase<double>);

private:
  LinearProblem MpProb_;
  LinearProblem ApProb_;
  LinearProblem FpProb_;
  LinearSolver<double> MpSolver_;
  LinearSolver<double> ApSolver_;
  LinearSolver<double> FSolver_;
};

}


#endif
