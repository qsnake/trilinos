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

#include <stdlib.h>

#include <stdexcept>

#include "DenseLinAlgPack_random_vector.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"

void DenseLinAlgPack::seed_random_vector_generator( unsigned int s )
{
  srand(s);
}

void DenseLinAlgPack::random_vector( value_type l, value_type u, DVectorSlice* v )
{
  if(!v)
    throw std::invalid_argument( "random_vector(...) : Error, "
      "v can not be NULL" );
  if( l > u )
    throw std::invalid_argument( "random_vector(...) : Error, "
      "l can not be greater than u" );
  for( DVectorSlice::iterator itr = v->begin(); itr != v->end(); )
    *itr++ = l + (double(rand())/RAND_MAX) * (u -l);
}
