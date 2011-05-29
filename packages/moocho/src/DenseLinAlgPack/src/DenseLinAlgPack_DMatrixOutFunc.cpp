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

#include "DenseLinAlgPack_DMatrixOutFunc.hpp"
#include "DenseLinAlgPack_DVectorOutFunc.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "Teuchos_FancyOStream.hpp"

std::ostream& DenseLinAlgPack::output(std::ostream& out_arg, const DMatrixSlice& gms
  , LinAlgPackIO::fmtflags extra_flags )
{
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::getFancyOStream(Teuchos::rcp(&out_arg,false));
  Teuchos::OSTab tab(out);
  int w = out->width(0) - 1; // get the set width (minus 1 since a space is inserted)

  if( !(extra_flags & LinAlgPackIO::ignore_dim_bit) ) {
    std::ios_base::fmtflags old = out->flags();
    *out
      << std::setw(0) << std::left << gms.rows() << ' ' << gms.cols()
      << std::endl;
    out->flags(old);
  }

  if( gms.rows() && gms.cols() ) {
    for(size_type i = 1; i <= gms.rows();++i) {
      const DVectorSlice& vs =gms.row(i); 
      DVectorSlice::const_iterator itr = vs.begin();
      for( size_type j = 1; itr != vs.end(); ++j, ++itr ) {
        *out << " " << std::setw(w) << (*itr) << ":" << i << ":" << j;
      }
      *out << std::endl;
    }
  }
  
  return out_arg;
}
