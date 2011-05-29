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

#ifndef SUNDANCE_NITSCHEBC_H
#define SUNDANCE_NITSCHEBC_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"


namespace Sundance
{
class CellFilter;
class QuadratureFamily;

/**
 * This function forms the expressions that apply the Dirichlet 
 * BC \f$u=u_{BC}\f$ via Nitsche's method for the Poisson operator.
 */
Expr NitschePoissonDirichletBC(int dim,
  const CellFilter& cells,
  const QuadratureFamily& quad,
  const Expr& kappa,
  const Expr& v,
  const Expr& u,
  const Expr& uBC,
  const double& gamma);

/**
 * This function forms the expressions that apply the Dirichlet 
 * BC \f${\bf u}={\bf u}_{BC}\f$ via Nitsche's method for the Stokes operator.
 * 
 * \param cells the surface on which the BC is to be applied
 * \param quad the quadrature rule to be used
 * \param nu the viscosity, which may be a function of velocity
 * \param v velocity test function
 * \param q pressure test function
 * \param u velocity unknown function
 * \param p pressure unknown function
 * \param uBC specified velocity at surface
 */
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
  );




}


#endif
