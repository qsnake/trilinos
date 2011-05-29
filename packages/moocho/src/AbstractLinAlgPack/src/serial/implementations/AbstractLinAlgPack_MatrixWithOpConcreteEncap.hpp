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

#ifndef MATRIX_WITH_OP_CONCRETE_ENCAP_H
#define MATRIX_WITH_OP_CONCRETE_ENCAP_H

#include "AbstractLinAlgPack_MatrixOp.hpp"

namespace AbstractLinAlgPack {

/** \brief This template class defines the storage for a concrete matrix
  * class that operations are based on.
  *
  * The default copy constructor and assignment operator are allowed.
  */
template<class M>
class MatrixWithOpConcreteEncap : public virtual MatrixOp
{
public:

  // /////////////////////////////////////////////////////
  /** @name Representation access */
  //@{

  /// The compiler did not generate this default constructor
  MatrixWithOpConcreteEncap()
  {}

  /// This constructor will have to be overridden.
  MatrixWithOpConcreteEncap(const M& m) : m_(m)
  {}

  /// Get the underlying M object
  M& m() {
    return m_;
  }

  /** \brief . */
  const M& m() const {
    return m_;
  }

  //@}	// end Representation access

  // /////////////////////////////////////////////////////
  // Overridden from Matrix

  /** \brief . */
  size_type rows() const;

  /** \brief . */
  size_type cols() const;

  // /////////////////////////////////////////////////////
  // Overridden from MatrixOp

  /** \brief . */
  MatrixOp& operator=(const MatrixOp& m);

private:
  M m_;

};	// end class MatrixWithOpConcreteEncap<M>

// Template definitions

template<class M>
size_type MatrixWithOpConcreteEncap<M>::rows() const {
  return m().rows();
}

template<class M>
size_type MatrixWithOpConcreteEncap<M>::cols() const {
  return m().cols();
}

template<class M>
MatrixOp& MatrixWithOpConcreteEncap<M>::operator=(const MatrixOp& m) {
  if(&m == this) return *this;	// assignment to self
  const MatrixWithOpConcreteEncap<M> *p_m = dynamic_cast<const MatrixWithOpConcreteEncap<M>*>(&m);
  if(p_m) {
    m_ = p_m->m_;
  }
  else {
    throw std::invalid_argument("MatrixWithOpConcreteEncap<M>::operator=(const MatrixOp& m)"
      " : The concrete type of m is not a subclass of MatrixWithOpConcreteEncap<M> as expected" );
  }
  return *this;
}

}	// end namespace AbstractLinAlgPack 

#endif	// MATRIX_WITH_OP_CONCRETE_ENCAP_H
