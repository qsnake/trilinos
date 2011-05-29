// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "ConstrainedOptPack_MeritFuncNLESqrResid.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"

namespace ConstrainedOptPack {

MeritFuncNLESqrResid::MeritFuncNLESqrResid()
  : deriv_(0.0)
{}

value_type MeritFuncNLESqrResid::calc_deriv( const Vector& c_k )
{
  using LinAlgOpPack::dot;
  return deriv_ = - dot(c_k,c_k);
}

// Overridden from MeritFuncNLP

value_type MeritFuncNLESqrResid::value(const Vector& c) const
{
  using LinAlgOpPack::dot;
  return 0.5 * dot(c,c);
}

value_type MeritFuncNLESqrResid::deriv() const
{
  return deriv_;
}

void MeritFuncNLESqrResid::print_merit_func(std::ostream& out
  , const std::string& L ) const
{
  out
    << L << "*** Define a square of constraint residuals merit funciton\n"
    << L << "*** (assumes Gc_k'*d_k + c_k = 0):\n"
    << L << "phi(c) = 1/2 * dot(c,c)\n"
    << L << "Dphi(x_k,d_k) = - dot(c_k,c_k)\n";
}

}	// end namespace ConstrainedOptPack 
