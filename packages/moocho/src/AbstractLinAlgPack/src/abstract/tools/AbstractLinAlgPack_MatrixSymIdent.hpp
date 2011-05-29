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

#ifndef ALAP_MATRIX_SYM_IDENTITY_H
#define ALAP_MATRIX_SYM_IDENTITY_H

#include "AbstractLinAlgPack_MatrixSymOpNonsing.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"

namespace AbstractLinAlgPack {

/** \brief Matrix subclass for a scaled identity matrix.
 *
 * More operations will be overridden as they are needed by various applications.
 */
class MatrixSymIdent : virtual public MatrixSymOpNonsing {
public:
  
  /** @name Constructors/initializers */
  //@{

  /// Calls <tt>this->initialize()</tt>.
  MatrixSymIdent(
    const VectorSpace::space_ptr_t&          vec_space = Teuchos::null
    ,const value_type                        scale     = 1.0
    );

  /** \brief . */
  void initialize(
    const VectorSpace::space_ptr_t&          vec_space
    ,const value_type                        scale       = 1.0
  );

  //@}

  /** @name Access */
  //@{

  /** \brief . */
  value_type scale() const;

  //@}

  /** @name Overridden from MatrixBase */
  //@{

  /// Returns 0 if not initalized.
  size_type rows() const;
  /// Returns <tt>this->rows()</tt>
  size_type nz() const;

  //@}

  /** @name Overridden from MatrixOp */
  //@{

  /** \brief . */
  const VectorSpace& space_cols() const;
  /** \brief . */
  std::ostream& output(std::ostream& out) const;
  /** \brief . */
  void Vp_StMtV(
    VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    ,const Vector& v_rhs2, value_type beta ) const;

  //@}

  /** @name Overridden from MatrixNonsing */
  //@{

  /** \brief . */
  void V_InvMtV(
    VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
    ,const Vector& v_rhs2 ) const;

  //@}

private:

  VectorSpace::space_ptr_t  vec_space_;
  value_type                scale_;
  
}; // end class MatrixSymIdent

// ///////////////////////////////////////////
// Inline members

inline
value_type MatrixSymIdent::scale() const
{
  return scale_;
}

} // end namespace AbstractLinAlgPack

#endif // ALAP_MATRIX_SYM_IDENTITY_H
