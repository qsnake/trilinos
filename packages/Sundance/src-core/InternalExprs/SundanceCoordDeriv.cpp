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

#include "SundanceCoordDeriv.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceDeriv.hpp"


using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;

CoordDeriv::CoordDeriv(int dir)
  : DerivBase(), dir_(dir)
{}

bool CoordDeriv::lessThan(const Deriv& other) const
{
  /* First compare types: spatial derivs are ranked lower than 
   * functional derivs, so if the other guy is a func deriv, I'm lower */
  if (other.isFunctionalDeriv()) return true;

  /* we can now safely cast the other guy to a spatial deriv and compare
   * directions */
  const CoordDeriv* s = other.coordDeriv();

  return dir_ < s->dir_;
}

string CoordDeriv::toString() const 
{
  return CoordExpr::coordName(dir_, "");
}


