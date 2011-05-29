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
#include <iomanip>

#include "DenseLinAlgPack_DVectorOutFunc.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"

std::ostream& DenseLinAlgPack::output(std::ostream& os, const DVectorSlice& vs
  , LinAlgPackIO::fmtflags extra_flags)
{
  int w = os.width(0) - 1; // get the set width (minus 1 since a space is inserted)

  if( !(extra_flags & LinAlgPackIO::ignore_dim_bit) )
    os << std::setw(0) << std::left << vs.dim() << std::endl << std::right;

  DVectorSlice::const_iterator itr = vs.begin();
  for( size_type i = 1; itr != vs.end(); ++i, ++itr ) {
    os << " " << std::setw(w) << (*itr) << ":" << i; // insert a space to be sure there is white space
                                                     // inbetween adjacent elements.
  }

  if( !(extra_flags & LinAlgPackIO::no_insert_newlines_bit) )
    os << std::endl;

  return os;
}
