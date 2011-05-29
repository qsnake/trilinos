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

#include "DenseLinAlgPack_DVectorInFunc.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"

namespace {	// Local implementation
std::istream& input_vs(std::istream& is, DenseLinAlgPack::DVectorSlice* vs, const char func[]);
}

std::istream& DenseLinAlgPack::input(std::istream& is, DVector* v, LinAlgPackIO::fmtflags extra_flags) {
  if( !(extra_flags & LinAlgPackIO::ignore_dim_bit) ) {
    size_type n;
    is >> n;
    if(is.fail())
      throw LinAlgPackIO::InputException("DenseLinAlgPack::input() {DVector}:  Input operation of vector dimension failed.  Check that the constant n is a valid integer.");
    if(is.bad())
      throw std::ios_base::failure("DenseLinAlgPack::input() {DVector}: Input operation failed because the stream became currupted.");
    v->resize(n);
  }
  return input_vs(is,&(*v)(),"DenseLinAlgPack::input() {DVector}");
}

std::istream& DenseLinAlgPack::input(std::istream& is, DVectorSlice* vs, LinAlgPackIO::fmtflags extra_flags) {
  if( !(extra_flags & LinAlgPackIO::ignore_dim_bit) ) {
    size_type n;
    is >> n;
    if(is.fail())
      throw LinAlgPackIO::InputException("DenseLinAlgPack::input() {DVectorSlice}:  Input operation of vector dimension failed.  Check that the constant n is a valid integer.");
    if(is.bad())
      throw std::ios_base::failure("DenseLinAlgPack::input() {DVectorSlice}: Input operation failed because the stream became currupted.");
    DenseLinAlgPack::Vp_V_assert_sizes( vs->dim(), n );
  }
  return input_vs(is,vs,"DenseLinAlgPack::input() {DVectorSlice}");
}


// ///////////////////////////////////
// Local implementation

namespace {

// Read in a specified number of elements into a DVectorSlice object.
// The dim of vs is not checked.  If an element input operation fails or the end of the file
// is reached before all of the elements are read in then a LinAlgPackIO::InputException is thrown.
// If the stream becomes currupted durring the input then a std::ios_base::failure exception
// is thrown.  The state of the input steam remains the same on return accept for the char's
// that have been extracted.
std::istream& input_vs(std::istream& is, DenseLinAlgPack::DVectorSlice* vs, const char func[]) {
  using std::ios_base;
  using DenseLinAlgPack::DVectorSlice;
  if(!vs->dim()) return is;	// If there are no elements to read in just return
  ios_base::iostate old_state = is.exceptions();		// save the old state
  is.exceptions(ios_base::badbit | ios_base::failbit);
  try {
    // Read in the elements
    for(DVectorSlice::iterator itr = vs->begin(); itr != vs->end(); ++itr)
      is >> *itr;
  }
  catch(std::ios_base::failure& excpt) {
    is.exceptions(old_state);
    if(is.bad()) throw;	// The stream was bad so rethrow the exception
    if(is.fail()) {
      std::ostringstream os;
      os << func << ":  An vector element input failed.  Check that the vector element is a valid C number.  "
         << excpt.what();
      throw DenseLinAlgPack::LinAlgPackIO::InputException(os.str());			
    }
    if(is.eof()) {
      std::ostringstream os;
      os << func << ":  DVector input failed.  The end of the file was found before all of the elements where read in.  "
         << excpt.what();;
      throw DenseLinAlgPack::LinAlgPackIO::InputException(os.str());			
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
