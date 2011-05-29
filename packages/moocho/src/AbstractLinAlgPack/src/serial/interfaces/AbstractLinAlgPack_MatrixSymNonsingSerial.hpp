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

#ifndef SLAP_MATRIX_SYM_NONSINGULAR_SERIAL_H
#define SLAP_MATRIX_SYM_NONSINGULAR_SERIAL_H

#include "AbstractLinAlgPack_MatrixNonsingSerial.hpp"
#include "AbstractLinAlgPack_MatrixSymNonsing.hpp"

namespace AbstractLinAlgPack {

/** \brief Abstract base class for all serial polymorphic symmetrix nonsingular matrices that
 * can be used to solve for linear systems relatively efficiently.
 *
 * The methods of this interface should not be called directly but instead through
 * the \ref MatrixSymNonsingularSerial_funcs "provided nonmember functions".
 */
class MatrixSymNonsingSerial
  : virtual public MatrixNonsingSerial
  , virtual public AbstractLinAlgPack::MatrixSymNonsing // doxygen needs full name
{
public:

  /** \brief . */
  using MatrixSymNonsing::M_StMtInvMtM;

  /** @name Level-3 */
  //@{

  /** \brief sym_gms_lhs = alpha * op(mwo) * inv(M) * op(mwo)'.
    *
    * The default implementation is based on the operation M_StInvMtM(...)
    * assuming that this \c M is a symmetric matrix.  For an efficient implementation
    * (for this = L*L' for instance) the subclass may want to override this function.
    */
  virtual void M_StMtInvMtM(
    DMatrixSliceSym* sym_gms_lhs, value_type alpha
    ,const MatrixOpSerial& mwo, BLAS_Cpp::Transp mwo_trans
    ,EMatrixDummyArg
    ) const;

  //@}

  /** @name Overridden from MatrixSymNonsing */
  //@{

  void M_StMtInvMtM(
    MatrixSymOp* sym_lhs, value_type alpha
    ,const MatrixOp& mwo, BLAS_Cpp::Transp mwo_trans
    ,EMatrixDummyArg
    ) const;

  //@}

};	// end class MatrixSymNonsingSerial

/** \defgroup MatrixSymNonsingularSerial_funcs MatrixSymNonsingSerial nonmember inline functions.
 *
 * These nonmember functions allow operations to be called on \c MatrixSymNonsingSerial objects
 * in similar manner to those in \c DenseLinAlgPack.
 */
//@{

inline
/// sym_gms_lhs = alpha * op(mwo) * inv(mswof) * op(mwo)'
void M_StMtInvMtM(
  DMatrixSliceSym* sym_gms_lhs, value_type alpha
  ,const MatrixOpSerial& mwo, BLAS_Cpp::Transp mwo_trans
  ,const MatrixSymNonsingSerial& mswons
  ,MatrixSymNonsingSerial::EMatrixDummyArg mwo_rhs
   )
{
  mswons.M_StMtInvMtM(sym_gms_lhs,alpha,mwo,mwo_trans,mwo_rhs);
}

//@}

} // end namespace AbstractLinAlgPack

#endif	// SLAP_MATRIX_SYM_NONSINGULAR_SERIAL_H
