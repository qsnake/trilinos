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

#include "SundanceNitscheBC.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceVectorCalculus.hpp"
#include "SundanceCellVectorExpr.hpp"
#include "SundanceCellDiameterExpr.hpp"

namespace Sundance
{

Expr NitschePoissonDirichletBC(int dim,
  const CellFilter& cells,
  const QuadratureFamily& quad,
  const Expr& kappa,
  const Expr& v,
  const Expr& u,
  const Expr& uBC,
  const double& gamma)
{
  Expr grad = gradient(dim);
  Expr n = CellNormalExpr(dim,"n");
  Expr h = new CellDiameterExpr();

  return Integral(cells, -kappa*(v*((n*grad)*u) + (u-uBC)*((n*grad)*v))
    + gamma/h * v*(u-uBC), quad);
}

Expr NitscheStokesNoSlipBC(const CellFilter& cells,
  const QuadratureFamily& quad,
  const Expr& nu,
  const Expr& v,
  const Expr& q,
  const Expr& u,
  const Expr& p,
  const Expr& uBC,
  const double& gamma1,
  const double& gamma2
  )
{
  TEST_FOR_EXCEPT(nu.size() != 1);
  TEST_FOR_EXCEPT(v.size() != u.size());
  TEST_FOR_EXCEPT(uBC.size() != u.size());

  int dim = uBC.size();

  Expr grad = gradient(dim);
  Expr n = CellNormalExpr(dim,"n");
  Expr h = new CellDiameterExpr();
  
  return Integral(cells,
    -nu*(v*((n*grad)*(u-uBC)) + (u-uBC)*((n*grad)*v))
    + p*(v*n) + q*((u-uBC)*n)
    + gamma1/h * v*(u-uBC) + gamma2/h * (v*n) * ((u-uBC)*n),
    quad);
}



}
