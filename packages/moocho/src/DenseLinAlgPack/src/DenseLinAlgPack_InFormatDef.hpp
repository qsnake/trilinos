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
//
// Template definition file.

#ifndef LINALGPACK_IN_FORMAT_DEF_H
#define LINALGPACK_IN_FORMAT_DEF_H

#include <istream>

#include "DenseLinAlgPack_InFormatDecl.hpp"

namespace DenseLinAlgPack {

template<class T>
std::istream& operator>>(std::istream& is, const LinAlgPackIO::bound_format<T>& bf) {
  using LinAlgPackIO::ios_format_memento;

  ios_format_memento old_format = ios_format_memento::save_format(is);

  try {
    bf.f().set_format(is);
    input( is, &const_cast< LinAlgPackIO::bound_format<T>&>(bf).obj()
         , bf.f().extra_flags().flags() );
  }
  catch(...) {
    old_format.set_format(is);
    throw;
  }

  old_format.set_format(is);
  return is;
}

}	// end namespace DenseLinAlgPack

#endif // LINALGPACK_IN_FORMAT_DEF_H
