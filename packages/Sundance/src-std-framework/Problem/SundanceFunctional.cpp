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

#include "SundanceFunctional.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceTestFunction.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceEssentialBC.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceListExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceAssembler.hpp"


using namespace Sundance;
using namespace Teuchos;
using namespace TSFExtended;


Functional::Functional(const Mesh& mesh, const Expr& integral, 
  const TSFExtended::VectorType<double>& vecType)
  : mesh_(mesh),
    integral_(integral),
    bc_(),
    vecType_(vecType)
{
  
}

Functional::Functional(
  const Mesh& mesh, 
  const Expr& integral,  
  const Expr& essentialBC,
  const TSFExtended::VectorType<double>& vecType)
  : mesh_(mesh),
    integral_(integral),
    bc_(essentialBC),
    vecType_(vecType)
{
  
}


LinearProblem Functional::linearVariationalProb(
  const Expr& var,
  const Expr& varEvalPts,
  const Expr& unk,
  const Expr& fixed,
  const Expr& fixedEvalPts) const 
{

  Array<Expr> zero(unk.size());
  for (int i=0; i<unk.size(); i++) 
    {
      Expr z = new ZeroExpr();
      zero[i] = z;
    }

  Expr unkEvalPts = new ListExpr(zero);

  Expr unkParams;
  Expr fixedParams;
  Expr unkParamValues;
  Expr fixedParamValues;

  RCP<EquationSet> eqn 
    = rcp(new EquationSet(integral_, bc_, 
                          tuple(var), tuple(varEvalPts),
                          tuple(unk), tuple(unkEvalPts), 
                          unkParams, unkParamValues,
                          tuple(fixed), tuple(fixedEvalPts)));

  RCP<Assembler> assembler 
    = rcp(new Assembler(mesh_, eqn, tuple(vecType_), tuple(vecType_), false));

  return LinearProblem(assembler);
}

NonlinearProblem Functional
::nonlinearVariationalProb(const Expr& var,
                           const Expr& varEvalPts,
                           const Expr& unk,
                           const Expr& unkEvalPts,
                           const Expr& fixed,
                           const Expr& fixedEvalPts) const
{

  Expr unkParams;
  Expr fixedParams;
  Expr unkParamValues;
  Expr fixedParamValues;

  RCP<EquationSet> eqn 
    = rcp(new EquationSet(integral_, bc_, 
                          tuple(var), tuple(varEvalPts),
                          tuple(unk), tuple(unkEvalPts), 
                          fixedParams, fixedParamValues,
                          tuple(fixed), tuple(fixedEvalPts)));

  RCP<Assembler> assembler 
    = rcp(new Assembler(mesh_, eqn, tuple(vecType_), tuple(vecType_), false));

  return NonlinearProblem(assembler, unkEvalPts);
}

FunctionalEvaluator Functional::evaluator(const Expr& var,
                                          const Expr& varEvalPts,
                                          const Expr& fixed,
                                          const Expr& fixedEvalPts) const 
{

  Expr unkParams;
  Expr fixedParams;
  Expr unkParamValues;
  Expr fixedParamValues;

  return FunctionalEvaluator(mesh_, integral_, bc_,
                             var, 
                             varEvalPts, 
                             fixed, fixedEvalPts,
                             vecType_);
}


FunctionalEvaluator Functional::evaluator(const Expr& var,
                                          const Expr& varEvalPts) const 
{
  return FunctionalEvaluator(mesh_, integral_, bc_,
                             var, 
                             varEvalPts,
                             vecType_);
}


namespace Sundance
{

double L2Norm(const Mesh& mesh, const CellFilter& domain,
  const Expr& f, const QuadratureFamily& quad)
{
  Expr I2 = Integral(domain, f*f, quad);

  return sqrt(evaluateIntegral(mesh, I2));
}


double H1Seminorm(
  const Mesh& mesh,
  const CellFilter& filter,
  const Expr& f,
  const QuadratureFamily& quad)
{
  Expr grad = gradient(mesh.spatialDim());
  return L2Norm(mesh, filter, grad*f, quad);
}
  
double H1Norm(
  const Mesh& mesh,
  const CellFilter& filter,
  const Expr& f,
  const QuadratureFamily& quad)
{
  Expr grad = gradient(mesh.spatialDim());
  Expr g = grad*f;
  Expr I2 = Integral(filter, f*f + g*g, quad);

  return sqrt(evaluateIntegral(mesh, I2));  
}
}
