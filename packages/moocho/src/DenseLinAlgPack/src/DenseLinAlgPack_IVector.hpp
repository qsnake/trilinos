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

#ifndef IVECTOR_H
#define IVECTOR_H

#include <assert.h>

#include <valarray>

#include "DenseLinAlgPack_Types.hpp"
#include "Teuchos_TestForException.hpp"

namespace DenseLinAlgPack {
/** \brief . */
/* * Fortran compatable integer vector for holding the  pivot information for
 * the elements of a vector, or the rows or columns of a matrix.
 */
class IVector : public std::valarray<DenseLinAlgPack::size_type> {
public:

  // STL typedefs
  typedef DenseLinAlgPack::index_type		value_type;
  typedef DenseLinAlgPack::size_type		size_type;
  typedef value_type&					reference;
  typedef const value_type&			const_reference;
  typedef value_type*					iterator;
  typedef const value_type*			const_iterator;
  typedef std::valarray<size_type>	valarray;

  // constructors

  /** \brief . */
  IVector();
  /** \brief . */
  IVector(size_type n);
  /** \brief . */
  IVector(const value_type& val, size_type n);
  /** \brief . */
  IVector(const value_type* p, size_type n);

  /// Resize on assignment
  IVector& operator=(const IVector&);

  /// 1-based element access (range checked if TEUCHOS_DEBUG is defined)
  reference operator()(size_type i);
  /// 1-based element access (range checked if TEUCHOS_DEBUG is defined)
  const_reference operator()(size_type i) const;

  /// STL iterator
  iterator begin();
  /// STL iterator
  const_iterator begin() const;
  /// STL iterator
  iterator end();
  /// STL iterator
  const_iterator end() const;

}; // end class IVector

// Inline definitions

inline IVector::IVector() : std::valarray<size_type>()
{}

inline IVector::IVector(size_type n) : std::valarray<size_type>(n)
{}

inline IVector::IVector(const value_type& val, size_type n) : std::valarray<size_type>(val,n)
{}

inline IVector::IVector(const value_type* p, size_type n) : std::valarray<size_type>(p,n)
{}

inline IVector& IVector::operator=(const IVector& iv)
{
  this->resize(iv.size());
  std::valarray<DenseLinAlgPack::size_type>::operator=(iv);
  return *this;
}

inline IVector::reference IVector::operator()(size_type i)
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !(  1 <= i && i <= static_cast<size_type>(size())  ) );
#endif
  return operator[](i-1);
}

inline IVector::const_reference IVector::operator()(size_type i) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !(  1 <= i && i <= static_cast<size_type>(size())  ) );
#endif
  return const_cast<IVector*>(this)->operator[](i-1);
}

inline IVector::iterator IVector::begin()
{	return &operator[](0); }

inline IVector::const_iterator IVector::begin() const
{	return &(const_cast<IVector*>(this)->operator[](0)); }

inline IVector::iterator IVector::end()
{	return begin() + size(); }

inline IVector::const_iterator IVector::end() const
{	return begin() + size(); }

}	// end namespace DenseLinAlgPack

#endif // IVECTOR_H
