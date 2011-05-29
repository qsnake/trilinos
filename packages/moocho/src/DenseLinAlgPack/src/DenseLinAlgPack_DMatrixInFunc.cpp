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

#include <sstream>

#include "InputStreamHelperPack_EatInputComment.hpp"
#include "DenseLinAlgPack_DMatrixInFunc.hpp"
#include "DenseLinAlgPack_DVectorInFunc.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"

namespace {	// Local inplementation
std::istream& input_gms(std::istream& is, DenseLinAlgPack::DMatrixSlice* gms, const char func[]);
}

std::istream& DenseLinAlgPack::input(std::istream& is, DMatrix* gm, LinAlgPackIO::fmtflags extra_flags) {
  if( !(extra_flags & LinAlgPackIO::ignore_dim_bit) ) {
    size_type m, n;
    is >> m >> n;
    if(is.fail())
      throw LinAlgPackIO::InputException( "DenseLinAlgPack::input() {DMatrix}: "
        "Input operation of matrix dimension failed.  Check that the constant n "
        "is a valid integer." );
    if(is.bad())
      throw std::ios_base::failure( "DenseLinAlgPack::input() {DMatrix}: "
        "Input operation failed because the stream became currupted." );
    gm->resize(m,n);
  }
  DMatrixSlice gms = (*gm)();
  return input_gms(is,&gms,"DenseLinAlgPack::input() {DMatrix}");
}

std::istream& DenseLinAlgPack::input(std::istream& is, DMatrixSlice* gms, LinAlgPackIO::fmtflags extra_flags) {
  if( !(extra_flags & LinAlgPackIO::ignore_dim_bit) ) {
    size_type m, n;
    is >> m >> n;
    if(is.fail())
      throw LinAlgPackIO::InputException( "DenseLinAlgPack::input() {DMatrixSlice}: "
        "Input operation of matrix dimension failed.  Check that the constant n "
        " is a valid integer.");
    if(is.bad())
      throw std::ios_base::failure( "DenseLinAlgPack::input() {DMatrixSlice}: "
        "Input operation failed because the stream became currupted." );
    DenseLinAlgPack::assert_gms_lhs(*gms,m,n);
  }
  return input_gms( is, gms, "DenseLinAlgPack::input() {DMatrixSlice}" );
}

// //////////////////////
// Local implementation

namespace {

// Read in a specified number of elements into a DMatrixSlice object.
// The dim of gms is not checked.  If an element input operation fails or the end of the file
// is reached before all of the elements are read in then a LinAlgPackIO::InputException is thrown.
// If the stream becomes currupted durring the input then a std::ios_base::failure exception
// is thrown.  The state of the input steam remains the same on return accept for the char's
// that have been extracted.
std::istream& input_gms(std::istream& is, DenseLinAlgPack::DMatrixSlice* gms, const char func[]) {
  using std::ios_base;
  using DenseLinAlgPack::size_type;
  using DenseLinAlgPack::DVectorSlice;
  if(!gms->rows()) return is;	// If we are inputting an unsized matrix then there are no elements
                // to extract so just return.
  ios_base::iostate old_state = is.exceptions();		// save the old state
  is.exceptions(ios_base::badbit | ios_base::failbit | ios_base::eofbit);
  try {
    // Read in the rows
    for(size_type i = 1; i <= gms->rows(); ++i) {
      InputStreamHelperPack::eat_comment_lines(is,'*');
      DVectorSlice gms_row_i = gms->row(i);
      DenseLinAlgPack::input( is, &gms_row_i
        , (DenseLinAlgPack::LinAlgPackIO::fmtflags)(DenseLinAlgPack::LinAlgPackIO::ignore_dim_bit) );
    }
  }
  catch(...) {
    is.exceptions(old_state);
    throw;
  }
  is.exceptions(old_state);
  return is;
}

}	// end namespace
