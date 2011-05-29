#include "SundancePCDPreconditioner.hpp"
#include "Sundance.hpp"
#include "TSFLinearSolverBuilder.hpp"
#include "TSFLinearSolverImpl.hpp"
#include "TSFSimpleIdentityOpDecl.hpp"
#include "TSFSimpleComposedOpDecl.hpp"
#include "TSFSimpleBlockOpDecl.hpp"
#include "TSFSimpleScaledOpDecl.hpp"
#include "TSFGenericRightPreconditioner.hpp"
#include "TSFLinearCombinationImpl.hpp"

using namespace TSFExtended;
using namespace Sundance;

PCDPreconditionerFactory::PCDPreconditionerFactory(
  const ParameterList& params,
  const LinearProblem& MpProb,
  const LinearProblem& ApProb,
  const LinearProblem& FpProb
  )
  : MpProb_(MpProb),
    ApProb_(ApProb),
    FpProb_(FpProb),
    MpSolver_(),
    ApSolver_(),
    FSolver_()
{
  ParameterList msParams = params.sublist("MpSolver");
  MpSolver_ = LinearSolverBuilder::createSolver(msParams);
  ParameterList asParams = params.sublist("ApSolver");
  ApSolver_ = LinearSolverBuilder::createSolver(asParams);
  ParameterList fsParams = params.sublist("FSolver");
  FSolver_ = LinearSolverBuilder::createSolver(fsParams);
}

Preconditioner<double> 
PCDPreconditionerFactory::
createPreconditioner(const LinearOperator<double>& K) const
{
  Tabs tab;

  LinearOperator<double> F = K.getBlock(0,0);
  F.setName("F");
  LinearOperator<double> FInv = inverse(F, FSolver_);
  FInv.setName("FInv");
  LinearOperator<double> Bt = K.getBlock(0,1);
  Bt.setName("Bt");


  LinearOperator<double> Fp = FpProb_.getOperator();

  LinearOperator<double> Mp = MpProb_.getOperator();
  Mp.setName("Mp");

  LinearOperator<double> MpInv = inverse(Mp, MpSolver_);
  MpInv.setName("MpInv");

  LinearOperator<double> Ap = ApProb_.getOperator();
  Ap.setName("Ap");

  LinearOperator<double> ApInv = inverse(Ap, ApSolver_);
  ApInv.setName("ApInv");


  VectorSpace<double> pDomain = Bt.domain();
  VectorSpace<double> uDomain = F.domain();

  LinearOperator<double> Iu = identityOperator(uDomain);
  Iu.setName("Iu");
  LinearOperator<double> Ip = identityOperator(pDomain);
  Ip.setName("Ip");

  LinearOperator<double> XInv = MpInv * Fp * ApInv;

  VectorSpace<double> rowSpace = K.range();
  VectorSpace<double> colSpace = K.domain();
   
  LinearOperator<double> Q1 = makeBlockOperator(colSpace, rowSpace);
  Q1.setName("Q1");
  LinearOperator<double> Q2 = makeBlockOperator(colSpace, rowSpace);
  Q2.setName("Q2");
  LinearOperator<double> Q3 = makeBlockOperator(colSpace, rowSpace);
  Q3.setName("Q3");
   
  Q1.setBlock(0, 0, FInv);
  Q1.setBlock(1, 1, Ip);
  Q1.endBlockFill();
   
  Q2.setBlock(0, 0, Iu);
  Q2.setBlock(0, 1, -1.0*Bt);
  Q2.setBlock(1, 1, Ip);
  Q2.endBlockFill();
   
  Q3.setBlock(0, 0, Iu);
  Q3.setBlock(1, 1, -1.0*XInv);
  Q3.endBlockFill();
   
  LinearOperator<double> P1 = Q2 * Q3;
  LinearOperator<double> PInv = Q1 * P1;

  return new GenericRightPreconditioner<double>(PInv);
}


