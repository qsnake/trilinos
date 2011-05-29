#if 0

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

#include "AbstractLinAlgPack_SparseCOOReadMatrix.hpp"

// Throw an exception if the char is not ':'
namespace {
inline void assert_sep_char(char c) {
  if(c != ':')
    throw AbstractLinAlgPack::InputException("Sparse COO matrix input stream error:  The seperator between the element, row indice and column indice must be a \':\'");
}
inline void assert_eof(std::istream& istrm) {
  if(istrm.eof())
    throw AbstractLinAlgPack::InputException("Sparse COO matrix input stream error:  Premature end to the input file.");
}
}

void AbstractLinAlgPack::read_coo_into_valarrays(std::istream& istrm, size_type& m, size_type& n, size_type& nz
  , std::valarray<value_type>& a, std::valarray<indice_type>& ivect
  , std::valarray<indice_type>& jvect)
{
  // read in dimensions and resize
  istrm >> m;		assert_eof(istrm);
  istrm >> n;		assert_eof(istrm);
  istrm >> nz;
  a.resize(nz);
  ivect.resize(nz);
  jvect.resize(nz);

  // Read in the non-zero elements
  value_type	*p_a =			&a[0],
        *p_a_last =		p_a + nz;
  indice_type	*p_ivect =		&ivect[0],
        *p_jvect =		&jvect[0];

  for(; p_a != p_a_last; ++p_a, ++p_ivect, ++p_jvect) {
    const int bs = 50;
    char num[bs];
    char c;
    assert_eof(istrm);
    istrm.get(num, bs-1, ':');  assert_eof(istrm);  *p_a = ::atof(num);	// Read in ak 
    istrm.get(c); assert_eof(istrm); assert_sep_char(c);				// Read in ':'
    istrm.get(num, bs-1, ':'); assert_eof(istrm); *p_ivect = ::atoi(num);// Read in ik 
    istrm.get(c); assert_eof(istrm); assert_sep_char(c);				// Read in ':'
    istrm >> *p_jvect;										// Read in jk
  }
}

#endif // 0
