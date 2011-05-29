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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-13079
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "RTOpPack_print_sub_vector.hpp"
#include "Teuchos_FancyOStream.hpp"

std::ostream& RTOpPack::output(
  std::ostream& o_arg, const SubVector& v
  ,bool print_dim , bool newline
  )
{
  Teuchos::RCP<Teuchos::FancyOStream> o = Teuchos::getFancyOStream(Teuchos::rcp(&o_arg,false));
  //Teuchos::OSTab tab(o);
  int w = o->width(0) - 1; // get the set width (minus 1 since a space is inserted)
  if( print_dim )
    *o << std::setw(0) << std::left << v.subDim() << std::endl << std::right;
  //*o << std::setiosflags(std::ios::left) << std::setw(0) << v.subDim() 
  //	  << std::endl << std::setiosflags(std::ios::right);
  // RAB: 20030916: ToDo: Fix the above by hacking std::left and std::right in config header!
  const RTOp_value_type  *v_val        = v.values();
  const ptrdiff_t        v_val_s       = v.stride();
  for( RTOp_index_type i = 1; i <= v.subDim(); ++i, v_val+=v_val_s ) {
    // insert a space to be sure there is white space
    // inbetween adjacent elements.
    *o << " " << std::setw(w) << (*v_val) << ":" << i + v.globalOffset();
  }
  if(newline) *o << std::endl;
  return o_arg;
}
