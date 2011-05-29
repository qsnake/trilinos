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

#include "AbstractLinAlgPack_GenPermMatrixSliceOut.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSlice.hpp"

std::ostream& AbstractLinAlgPack::operator<<(
  std::ostream& out, const GenPermMatrixSlice& P
  )
{
  out	<< P.rows() << " " << P.cols() << " " << P.nz() << std::endl;
  if( P.is_identity() ) {
    out << "identity matrix\n";
  }
  else if( P.nz() ) {
    for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr )
      out << " " << itr->row_i() << ":" << itr->col_j();
    out << std::endl;
  }
  return out;
}
