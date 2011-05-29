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

#ifndef MATRIX_IDENT_CONCAT_H
#define MATRIX_IDENT_CONCAT_H

#include "ConstrainedOptPack_Types.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"

namespace ConstrainedOptPack {

/** \brief Matrix class for a matrix vertically concatonated with an identity matrix {abstract}.
 *
 * Represents an interface for a matrix that represents:
 \verbatim

 M = [ alpha*op(D) ]
     [      I      ]
 where:
     D_rng = [1,rows(op(D))]
     I_rng = [rows(op(D))+1,rows(op(D))+cols(op(D))]

 or

 M = [      I      ]
     [ alpha*op(D) ]
 where:
     D_rng = [cols(op(D))+1,rows(op(D))+cols(op(D))]
     I_rng = [1,cols(op(D))]
 \endverbatim
 * and \c I is a <tt>op(D).cols() x op(D).cols()</tt> indentity matrix and
 * the full matrix \c M is of order <tt>(op(D).rows() + op(D).cols()) x op(D).cols()</tt>.
 */
class MatrixIdentConcat
  : virtual public AbstractLinAlgPack::MatrixOp
{
public:

  /** @name Access to representation.
   */
  //@{
  /** \brief . */
  virtual Range1D D_rng() const = 0;
  /** \brief . */
  virtual Range1D I_rng() const = 0;
  /** \brief . */
  virtual value_type alpha() const = 0;
  /** \brief . */
  virtual const MatrixOp& D() const = 0;
  /** \brief . */
  virtual BLAS_Cpp::Transp D_trans() const = 0;
  //@}

  /** @name Overridden from MatrixBase */
  //@{
  /** \brief . */
  size_type rows() const;
  /** \brief . */
  size_type cols() const;
  /** \brief . */
  size_type nz() const;
  //@}

  /** @name Overridden from MatrixOp */
  //@{
  /** \brief . */
  std::ostream& output(std::ostream& out) const;
  /** \brief . */
  void Vp_StMtV(
    VectorMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    ,const Vector& vs_rhs2, value_type beta
    ) const;
  /** \brief . */
  void Vp_StMtV(
    VectorMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    ,const SpVectorSlice& sv_rhs2, value_type beta
    ) const;
  //@}

}; // end class MatrixIdentConcat

} // end namespace ConstrainedOptPack

#endif // MATRIX_IDENT_CONCAT_H
