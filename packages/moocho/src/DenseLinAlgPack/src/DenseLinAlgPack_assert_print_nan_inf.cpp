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

#include <ostream>
#include <sstream>
#include <iomanip>

#include "DenseLinAlgPack_assert_print_nan_inf.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "check_nan_inf.h"

bool DenseLinAlgPack::assert_print_nan_inf( const value_type& val, char name[]
  , bool throw_excpt, std::ostream* out )
{
  
  if( RTOp_is_nan_inf(val) ) {
    std::ostringstream omsg;
    omsg
      << "The scalar \"" << name
      << "\" = " << val << " is not a valid bounded number";
    if(out)
      *out << omsg.str() << std::endl;
    if( throw_excpt ) {
      if(out)
        out->flush();	
      throw NaNInfException( "assert_print_nan_inf(...) : Error, "
        + omsg.str() );
    }
    return false;
  }
  return true;
}

bool DenseLinAlgPack::assert_print_nan_inf( const DVectorSlice& v, char name[]
  , bool throw_excpt, std::ostream* out )
{
  
  bool has_nan_or_inf = false;
  bool printed_header = false;

  for( DVectorSlice::const_iterator v_itr = v.begin(); v_itr != v.end(); ++v_itr ) {
    if( RTOp_is_nan_inf(*v_itr) ) {
      if(out) {
        if(!printed_header) {
          *out
            << "The vector \"" << name
            << "\" has the following NaN or Inf entries\n";
          printed_header = true;
        }
        *out
          << name << "(" << v_itr - v.begin() + 1 << ") = "
          << *v_itr << std::endl;
      }
      has_nan_or_inf = true;
    }
  }
  if( has_nan_or_inf && throw_excpt ) {
    if(out)
      out->flush();	
    std::ostringstream omsg;
    omsg
      << "assert_print_nan_inf(...) : Error, the vector named "
      << name << " has at least one element which is NaN or Inf";
    throw NaNInfException( omsg.str() );
  }

  return !has_nan_or_inf;
}

bool DenseLinAlgPack::assert_print_nan_inf( const DMatrixSlice& m, char name[]
  , bool throw_excpt, std::ostream* out )
{
  
  bool has_nan_or_inf = false;
  bool printed_header = false;

  for( size_type j = 1; j <= m.cols(); ++j ) {
    const DVectorSlice& v = m.col(j);
    for( DVectorSlice::const_iterator v_itr = v.begin(); v_itr != v.end(); ++v_itr ) {
      if( RTOp_is_nan_inf(*v_itr) ) {
        if(out) {
          if(!printed_header) {
            *out
              << "The matrix \"" << name
              << "\" has the following NaN or Inf entries\n";
            printed_header = true;
          }
          *out
            << name << "(" << v_itr - v.begin() + 1 << "," << j << ") = "
            << *v_itr << std::endl;
        }
        has_nan_or_inf = true;
      }
    }
  }
  
  if( has_nan_or_inf && throw_excpt ) {
    if(out)
      out->flush();	
    std::ostringstream omsg;
    omsg
      << "assert_print_nan_inf(...) : Error, the matrix named "
      << name << " has at least one element which is NaN or Inf";
    throw NaNInfException( omsg.str() );
  }

  return has_nan_or_inf;
}
