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

#ifndef COO_MATRIX_OUT_FUNC_H
#define COO_MATRIX_OUT_FUNC_H

#include <ostream>

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** \brief Output stream function for COOMatrix.
  *
  * The format for the output is:
  *
  * #m  n  nz#\\
  * #a1:i1:j1  a2:i2:j2 .... anz:inz:jnz#\\
  *
  * In the above format, each non-zero element is given as a three item pair:
  * value of the non-zero element, row indice (1-based) of the non-zero element,
  * and the column indice (1-based) of the non-zero element.  There are
  * no spaces between the numbers and the \':\' charachter and there is
  * one whitespace character between elements.  There is a new line character
  * added after all of the elements have been output.
  */
std::ostream& output(std::ostream& o, const COOMatrix& coom);

} // end namespace AbstractLinAlgPack

#endif // COO_MATRIX_OUT_FUNC_H
