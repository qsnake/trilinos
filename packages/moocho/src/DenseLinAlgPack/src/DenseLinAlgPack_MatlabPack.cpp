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

#include <limits>
#include <iomanip>
#include <ostream>

#include "DenseLinAlgPack_MatlabPack.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"

namespace DenseLinAlgPack {

std::ostream& MatlabPack::out( std::ostream& o, const char* name, const DVectorSlice& vs
  , BLAS_Cpp::Transp trans )
{
  int p = o.precision();
  o.precision( std::numeric_limits<value_type>::digits10 + 3 );
  try {
    o << name << " =  [ ";
    for( DVectorSlice::const_iterator itr = vs.begin(); itr != vs.end(); ++itr )
      o << ' ' << *itr << ";";
    o << "];" << (trans == BLAS_Cpp::no_trans ? ' ' : '\'' ) << std::endl;
  }
  catch(...) {
    o.precision(p);
    throw;
  }
  o.precision(p);
  return o;
}

std::ostream& MatlabPack::out( std::ostream& o, const char* name, const DMatrixSlice& gms
  , BLAS_Cpp::Transp trans )
{
  int p = o.precision();
  o.precision( std::numeric_limits<value_type>::digits10 + 3 );
  try {
    o << name << " =  [\n";
    for( size_type i = 1; i <= gms.rows(); ++i ) {
      const DVectorSlice vs = gms.row(i);
      for( DVectorSlice::const_iterator itr = vs.begin(); itr != vs.end(); ++itr )
        o << *itr << ", ";
      o << ";\n";
    }
    o << "];" << (trans == BLAS_Cpp::no_trans ? ' ' : '\'' ) << std::endl;
  }
  catch(...) {
    o.precision(p);
    throw;
  }
  o.precision(p);
  return o;
}

}	// end namespace DenseLinAlgPack
