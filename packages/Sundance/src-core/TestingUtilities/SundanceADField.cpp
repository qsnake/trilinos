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


#include "SundanceADField.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceExceptions.hpp"


using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace SundanceTesting;

using namespace Teuchos;
using namespace std;

ADField::ADField(const ADBasis& basis, const double& coeff)
  : basis_(basis), coeff_()
{
  double* x = new double;
  *x = coeff;
  coeff_ = rcp(x);
}

ADReal ADField::evaluate() const
{
  return *coeff_ * basis_.evaluate(evalPoint());
}

void ADField::setCoeff(const double& c)
{
  *coeff_ = c;
}

ADReal ADField::operator+(const ADReal& x) const
{
  return evaluate() + x;
}

ADReal ADField::operator+(const ADField& x) const
{
  return evaluate() + x.evaluate();
}

ADReal ADField::operator-(const ADField& x) const
{
  return evaluate() - x.evaluate();
}

ADReal ADField::operator*(const ADField& x) const
{
  return evaluate() * x.evaluate();
}


ADReal ADField::operator/(const ADField& x) const
{
  return evaluate() / x.evaluate();
}


ADReal ADField::operator-(const ADReal& x) const
{
  return evaluate() - x;
}

ADReal ADField::operator+(const double& x) const
{
  return evaluate() + x;
}


ADReal ADField::operator-(const double& x) const
{
  return evaluate() - x;
}

ADReal ADField::operator-() const
{
  return -evaluate();
}

ADReal ADField::operator*(const ADReal& x) const
{
  return evaluate() * x;
}


ADReal ADField::operator*(const double& x) const
{
  return evaluate() * x;
}



ADReal ADField::operator/(const ADReal& x) const
{
  return evaluate() / x;
}


ADReal ADField::operator/(const double& x) const
{
  return evaluate() / x;
}


ADReal ADField::reciprocal() const
{
  return 1.0/evaluate();
}
