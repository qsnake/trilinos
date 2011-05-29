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

#ifndef COP_MATRIX_SYM_IDENTITY_SERIAL_H
#define COP_MATRIX_SYM_IDENTITY_SERIAL_H

#include "ConstrainedOptPack_Types.hpp"
#include "AbstractLinAlgPack_MatrixExtractInvCholFactor.hpp"
#include "AbstractLinAlgPack_MatrixSymOpNonsingSerial.hpp"

namespace ConstrainedOptPack {

/** \brief Matrix class for a serial scaled identity matrix.
 *
 * More operations will be overridden as they are needed by various applications.
 */
class MatrixSymIdentitySerial
  : virtual public AbstractLinAlgPack::MatrixSymOpNonsingSerial // doxygen needs full name
  , virtual public AbstractLinAlgPack::MatrixExtractInvCholFactor
{
public:
  
    /** @name Constructors/initalizers */
  //@{

  /// Calls <tt>this->initalize()</tt>
  MatrixSymIdentitySerial( size_type size = 1, value_type scale = 1.0 );

  /** \brief . */
  void initialize( size_type size, value_type scale );

  //@}

  /** @name Access */
  //@{

  /** \brief . */
  value_type scale() const;

  //@}

  /** Overridden from MatrixBase */
  //@{

  /** \brief . */
  size_type rows() const;
  /** \brief . */
  size_type nz() const;

  //@}

  /** Overridden from MatrixOp */
  //@{

  /** \brief . */
  std::ostream& output(std::ostream& out) const;

  //@}

  /** Overridden from MatrixOpSerial */
  //@{

  /** \brief . */
  void Vp_StMtV(
    DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    ,const DVectorSlice& vs_rhs2, value_type beta) const;

  //@}

  /** @name Overridden from MatrixNonsingSerial */
  //@{

  /** \brief . */
  void V_InvMtV(
    DVectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1,const DVectorSlice& vs_rhs2 ) const;

  //@}

  /** @name Overridden from MatrixSymNonsingSerial */
  //@{

  /** \brief . */
  void M_StMtInvMtM(
    DMatrixSliceSym* sym_gms_lhs, value_type alpha
    ,const MatrixOpSerial& mwo, BLAS_Cpp::Transp mwo_trans
    ,EMatrixDummyArg
    ) const;

  //@}

  /** @name Overridden from MatrixExtractInvCholFactor */
  //@{

  /** \brief . */
  void extract_inv_chol( DMatrixSliceTriEle* InvChol ) const;

  //@}

private:

  size_type   size_;
  value_type  scale_;

  
}; // end class MatrixSymIdentitySerial

// //////////////////////////////////////
// Inline members

inline
value_type MatrixSymIdentitySerial::scale() const
{
  return scale_;
}

} // end namespace ConstrainedOptPack

#endif // COP_MATRIX_SYM_IDENTITY_SERIAL_H
