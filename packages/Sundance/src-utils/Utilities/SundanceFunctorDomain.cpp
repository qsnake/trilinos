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

#include "SundanceFunctorDomain.hpp"
#include "SundanceExceptions.hpp"

using namespace Sundance;
using namespace Teuchos;


FunctorDomain::FunctorDomain() {;}

double FunctorDomain::lowerBound() const
{
  TEST_FOR_EXCEPTION(true, InternalError, 
                     "FunctorDomain::lowerBound() called for a domain without "
                     "a lower bound");
  return 0.0;
}

double FunctorDomain::upperBound() const
{
  TEST_FOR_EXCEPTION(true, InternalError, 
                     "FunctorDomain::upperBound() called for a domain without "
                     "an upper bound");

  return 0.0;
}

double FunctorDomain::excludedPoint() const
{
  TEST_FOR_EXCEPTION(true, InternalError, 
                     "FunctorDomain::excludedPoint() called for a domain without "
                     "an excluded point");

  return 0.0;
}

UnboundedDomain::UnboundedDomain() {;}

PositiveDomain::PositiveDomain() {;}

BoundedDomain::BoundedDomain(const double& lower, const double& upper) 
  : FunctorDomain(), lower_(lower), upper_(upper) 
{;}

LowerBoundedDomain::LowerBoundedDomain(const double& lower)
  : FunctorDomain(), lower_(lower)
{;}

NonzeroDomain::NonzeroDomain() {;}

