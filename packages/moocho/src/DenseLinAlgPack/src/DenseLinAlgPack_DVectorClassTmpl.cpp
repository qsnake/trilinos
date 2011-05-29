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

#include <iomanip>

#include "DenseLinAlgPack_DVectorClassTmpl.hpp"
#include "Teuchos_TestForException.hpp"

#ifdef LINALGPACK_CHECK_SLICE_SETUP
DenseLinAlgPack::size_type DenseLinAlgPack::vector_validate_sized(size_type size)
{
  TEST_FOR_EXCEPTION(
    !size, std::invalid_argument
    ,"vector_validate_sized(...) : Error, A vector region can not be created from an unsized vector.\n"
    );
  return size;
}
#endif

#ifdef LINALGPACK_CHECK_RANGE
void DenseLinAlgPack::vector_validate_range(size_type ubound, size_type max_ubound)
{
  TEST_FOR_EXCEPTION(
    ubound > max_ubound, std::out_of_range
    ,"vector_validate_range(...) : The upper bound is out of range.\n");
}
#endif

#ifdef LINALGPACK_CHECK_RANGE
void DenseLinAlgPack::vector_validate_subscript(size_type size, size_type i)
{
  TEST_FOR_EXCEPTION(
    i < 1 || i > size, std::out_of_range
    ,"vector_validate_subscript(size,i) : Error, Subscript i out of bounds.\n");
}
#endif

#ifdef LINALGPACK_CHECK_RHS_SIZES
void DenseLinAlgPack::assert_vs_sizes(size_type size1, size_type size2)
{
  TEST_FOR_EXCEPTION(
    size1 != size2, std::length_error
    ,"assert_vs_sizes(...) : Error, size1 = " << size1 << " != size2 = " << size2 << "\n");
}
#endif
