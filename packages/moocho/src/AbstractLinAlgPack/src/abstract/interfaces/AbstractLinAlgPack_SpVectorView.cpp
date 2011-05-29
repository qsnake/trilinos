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

#include "AbstractLinAlgPack_SpVectorView.hpp"

namespace {

// Setup some template classes to check at complile time that
// the layout of SpVectorSlice::element_type is proper.
template<int N, class T>
class assert_compile_time {
public:
  assert_compile_time()
  {
    // This should not compile if instantiated with a type T that
    // is not an integer.  However, if the compiler checks this
    // function without instantiating it, it can not cause an error
    // because it does not know the type of T to see if the
    // conversion is legal or not.
    T d;
    static_cast<int*>(d);
  }
};

// Template specialization for no error
template<>
class assert_compile_time<0,double> {
public:
  assert_compile_time()
  {}
};

// Validate that there is an integer stride between values
assert_compile_time<
    ((int)sizeof(AbstractLinAlgPack::SpVectorSlice::element_type)
   % (int)sizeof(AbstractLinAlgPack::value_type))
  , double
  >
validate_value_stride;

// Validate that there is an integer stride between indexes
assert_compile_time<
    ((int)sizeof(AbstractLinAlgPack::SpVectorSlice::element_type)
   % (int)sizeof(AbstractLinAlgPack::index_type))
  , double
  >
validate_index_stride;

// Compute the stride between values
const int values_stride = (int)sizeof(AbstractLinAlgPack::SpVectorSlice::element_type)
  / (int)sizeof(AbstractLinAlgPack::value_type);

// Compute the stride between indexes
const int indices_stride = (int)sizeof(AbstractLinAlgPack::SpVectorSlice::element_type)
  / (int)sizeof(AbstractLinAlgPack::index_type);

} // end namespace

RTOpPack::SparseSubVector
AbstractLinAlgPack::sub_vec_view(
  const SpVectorSlice&   sv_in
  ,const Range1D&        rng_in
  )
{
  const Range1D        rng = RangePack::full_range(rng_in,1,sv_in.dim());
  const SpVectorSlice  sv = sv_in(rng);

  RTOpPack::SparseSubVector  sub_vec;

  if(!sv.nz()) {
    sub_vec.initialize(
      rng.lbound()-1  // global_offset
      ,rng.size()     // sub_dim
      ,0              // nz
      ,NULL           // vlaues
      ,1              // values_stride
      ,NULL           // indices
      ,1              // indices_stride
      ,0              // local_offset
      ,1              // is_sorted
      );
  }
  else {
    SpVectorSlice::const_iterator
      itr = sv.begin();
    TEST_FOR_EXCEPT( !( itr != sv.end() ) );
    const value_type  *values  = &itr->value();
    if( sv.dim() && sv.nz() == sv.dim() && sv.is_sorted() ) {
      sub_vec.initialize(
        rng.lbound()-1    // global_offset
        ,rng.size()       // sub_dim
        ,values           // values
        ,values_stride    // values_stride
        );
    }
    else {
      const value_type   *values  = &itr->value();
      const index_type   *indexes = &itr->index();
      sub_vec.initialize(
        rng.lbound()-1    // global_offset
        ,sv.dim()         // sub_dim
        ,sv.nz()          // sub_nz
        ,values           // values
        ,values_stride    // values_stride
        ,indexes          // indices
        ,indices_stride   // indices_stride
        ,sv.offset()      // local_offset
        ,sv.is_sorted()   // is_sorted
        );
    }
  }
  
  return sub_vec;
}
