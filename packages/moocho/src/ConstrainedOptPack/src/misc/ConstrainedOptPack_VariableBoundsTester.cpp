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

#include "ConstrainedOptPack_VariableBoundsTester.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"

namespace ConstrainedOptPack {

// public

VariableBoundsTester::VariableBoundsTester(
  value_type   warning_tol
  ,value_type  error_tol
  )
  :warning_tol_(warning_tol)
  ,error_tol_(error_tol)
{}

bool VariableBoundsTester::check_in_bounds(
  std::ostream* out, bool print_all_warnings, bool print_vectors
  ,const Vector& xL, const char xL_name[]
  ,const Vector& xU, const char xU_name[]
  ,const Vector& x,  const char x_name[]
  )
{
  using AbstractLinAlgPack::max_near_feas_step;

  if(out)
    *out
      << "\n*** Checking that variables are in bounds\n";

  VectorSpace::vec_mut_ptr_t zero = x.space().create_member(0.0);
  std::pair<value_type,value_type>
    u = max_near_feas_step( x, *zero, xL, xU, warning_tol() );
  if(u.first < 0.0) {
    if(out)
      *out << "\nWarning! the variables " << xL_name << " <= " << x_name << " <= " << xU_name
        << " are out of bounds by more than warning_tol = "	<< warning_tol() << "\n";
    u = max_near_feas_step( x, *zero, xL, xU, error_tol() );
    if(u.first < 0.0) {
      if(out)
        *out << "\nError! the variables " << xL_name << " <= " << x_name << " <= " << xU_name
          << " are out of bounds by more than error_tol = " << error_tol() << "\n";
      return false;
    }
  }
  return true;
}

}	// end namespace ConstrainedOptPack
