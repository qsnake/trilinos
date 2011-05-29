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

#include "ConstrainedOptPack_MeritFuncCalc1DQuadratic.hpp"
#include "ConstrainedOptPack_MeritFuncCalc.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "Teuchos_TestForException.hpp"

namespace ConstrainedOptPack {

MeritFuncCalc1DQuadratic::MeritFuncCalc1DQuadratic(
  const MeritFuncCalc&      phi
  ,size_type                p
  ,const_VectorWithOp_ptr   d[]
  ,VectorMutable*     x
  )
  : phi_(phi), p_(p), x_(x)
{
  TEST_FOR_EXCEPTION(
    !(1 <= p && p <= 2 ), std::invalid_argument
    ,"MeritFuncCalc1DQuadratic::MeritFuncCalc1DQuadratic(...) : Error! "
    "p = " << p << " must be in the range 1 <= p <= 2"
    );
  for( size_type i = 0; i <= p; ++i )
    d_[i] = d[i];
}

value_type MeritFuncCalc1DQuadratic::operator()(value_type alpha) const
{
  using AbstractLinAlgPack::Vp_StV;
  *x_ = *d_[0];
  value_type alpha_i = alpha;
  for( size_type i = 1; i <= p_; ++i, alpha_i *= alpha ) {
    Vp_StV( x_, alpha_i, *d_[i] );
  }
  return phi_( *x_ );
}

value_type  MeritFuncCalc1DQuadratic::deriv() const
{
  return phi_.deriv();
}

void MeritFuncCalc1DQuadratic::print_merit_func(
  std::ostream& out, const std::string& L
  ) const
{
  out	<< L << "*** MeritFuncCalc1DQuadratic\n"
    << L << "x = xo + alpha*d[1] + alpha^2*d[2]\n";
  phi_.print_merit_func( out, L );
}

}	// end namespace ConstrainedOptPack
