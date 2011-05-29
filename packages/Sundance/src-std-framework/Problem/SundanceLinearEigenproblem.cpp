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

#include "SundanceLinearEigenproblem.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceTestFunction.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceEssentialBC.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceListExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceGaussianQuadrature.hpp"

#include "TSFLinearCombinationDecl.hpp"
#include "TSFLinearCombinationImpl.hpp"
#include "TSFLinearOperatorDecl.hpp"
#include "TSFVectorDecl.hpp"
#include "TSFSimpleDiagonalOpDecl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorImpl.hpp"
#include "TSFLinearCombinationImpl.hpp"
#include "TSFLinearOperatorImpl.hpp"
#include "TSFSimpleDiagonalOpImpl.hpp"
#endif

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using namespace TSFExtended;
using namespace TSFExtendedOps;

static Time& normalizationTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("Eigenfunction normalization"); 
  return *rtn;
}

static Time& makeEigensystemTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("Building eigensystem stiffness matrix"); 
  return *rtn;
}

LinearEigenproblem::LinearEigenproblem(
  const Mesh& mesh, const Expr& eqn,
  const Expr& v, const Expr& u,
  const VectorType<double>& vecType)
  : lumpMass_(false),
    kProb_(),
    mProb_(),
    M_(),
    MUnlumped_(),
    discSpace_()
{
  Expr empty;
  
  kProb_ = LinearProblem(mesh, eqn, empty, v, u, vecType);
  discSpace_ = *(kProb_.solnSpace()[0]);
}    

LinearEigenproblem::LinearEigenproblem(
  const Mesh& mesh, const Expr& eqn,
  const Expr& v, const Expr& u,
  const VectorType<double>& vecType,
  bool lumpedMass)
  : lumpMass_(lumpedMass),
    kProb_(),
    mProb_(),
    M_(),
    MUnlumped_(),
    discSpace_()
{
  Expr empty;
  
  kProb_ = LinearProblem(mesh, eqn, empty, v, u, vecType);
  mProb_ = makeMassProb(mesh, empty, v, u, vecType);
  discSpace_ = *(kProb_.solnSpace()[0]);
  MUnlumped_ = mProb_.getOperator();
  if (lumpMass_)
  {
    M_ = lumpedOperator(MUnlumped_);
  }
  else
  {
    M_ = MUnlumped_;
  }
}    


LinearEigenproblem::LinearEigenproblem(
  const Mesh& mesh, const Expr& eqn,
  const Expr& massExpr,
  const Expr& v, const Expr& u,
  const VectorType<double>& vecType,
  bool lumpedMass)
  : lumpMass_(lumpedMass),
    kProb_(),
    mProb_(),
    M_(),
    MUnlumped_(),
    discSpace_()
{
  Expr bc;
  kProb_ = LinearProblem(mesh, eqn, bc, v, u, vecType);
  mProb_ = makeMassProb(mesh, massExpr, v, u, vecType);
  discSpace_ = *(kProb_.solnSpace()[0]);

  MUnlumped_ = mProb_.getOperator();
  if (lumpMass_)
  {
    M_ = lumpedOperator(MUnlumped_);
  }
  else
  {
    M_ = MUnlumped_;
  }
  
}    

LinearProblem LinearEigenproblem::makeMassProb(
  const Mesh& mesh,
  const Expr& massExpr,
  const Expr& v, const Expr& u,
  const VectorType<double>& vecType) const
{
  Expr eqn;
  
  CellFilter interior = new MaximalCellFilter();
  QuadratureFamily quad = new GaussianQuadrature( 4 );
  if (massExpr.ptr().get()==0)
  {
    eqn = Integral(interior, v*u, quad);
  }
  else
  {
    eqn = Integral(interior, massExpr, quad);
  }
  Expr bc;
  LinearProblem rtn(mesh, eqn, bc, v, u, vecType);
  return rtn;
}


Array<Expr> LinearEigenproblem::makeEigenfunctions(
  Array<Vector<double> >& ev) const 
{
  TimeMonitor timer(normalizationTimer());

  Array<Expr> x(ev.size());
  CellFilter interior = new MaximalCellFilter();
  QuadratureFamily q = new GaussianQuadrature(2);
  for (int i=0; i<ev.size(); i++) 
  {
    x[i] = new DiscreteFunction(discSpace_, ev[i], "ev[" + Teuchos::toString(i)+"]");
    double N = 1.0;
    if (MUnlumped_.ptr().get())
    {
      N = ev[i] * (MUnlumped_ * ev[i]);
    }
    else
    {
      N = evaluateIntegral(discSpace_.mesh(), 
        Integral(interior, x[i]*x[i], q));
    }
    ev[i].scale(1.0/sqrt(N));
  }

  return x;
}


LinearOperator<double> 
LinearEigenproblem::lumpedOperator(const LinearOperator<double>& M) const 
{
  Vector<double> ones = M.domain().createMember();
  ones.setToConstant(1.0);
  Vector<double> m = M * ones;
  LinearOperator<double> rtn = diagonalOperator(m);

  return rtn;
}


Eigensolution LinearEigenproblem::solve(const Eigensolver<double>& solver) const 
{
  Array<std::complex<double> > ew;
  Array<Vector<double> > ev;

  LinearOperator<double> K;
  {
    TimeMonitor timer(makeEigensystemTimer());
    K = kProb_.getOperator();
  }
  
  solver.solve(K, M_, ev, ew);

  Array<Expr> eigenfuncs = makeEigenfunctions(ev);

  return Eigensolution(eigenfuncs, ew);
}

