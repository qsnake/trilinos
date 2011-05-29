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

#ifndef SLAP_MATRIX_SYM_WITH_OP_SERIAL_H
#define SLAP_MATRIX_SYM_WITH_OP_SERIAL_H

#include "AbstractLinAlgPack_MatrixOpSerial.hpp"
#include "AbstractLinAlgPack_MatrixSymOp.hpp"

namespace AbstractLinAlgPack {

/** \brief Abstract base class for all <tt>AbstractLinAlgPack::MatrixSymOp</tt> objects
 * implemented in shared memory space.
 *
 * This base class does a mapping from fully abstract linear algebra to shared memory
 * linear algebra.
 *
 * These methods should not be called directly but instead should be called through
 * the line \ref MatrixSymWithOpSerial_funcs "non-member functions" that are provided.
 */
class MatrixSymOpSerial
  : virtual public MatrixOpSerial
  , virtual public AbstractLinAlgPack::MatrixSymOp // doxygen needs full name
{
public:

  /** \brief . */
  using MatrixSymOp::Mp_StPtMtP;

  /** \brief sym_lhs = alpha * op(gpms_rhs') * M * op(gpms_rhs) + beta * sym_lhs.
    *
    * The default operation is based on <tt>this->Vp_StMtV(...)</tt> and assumes
    * that the matrix is symmetric.  Of course, a more efficient implementation
    * is often needed and the sublcass would like to override this.
    */
  virtual void Mp_StPtMtP(
    DMatrixSliceSym* sym_lhs, value_type alpha
    ,EMatRhsPlaceHolder dummy_place_holder
    ,const GenPermMatrixSlice& gpms_rhs, BLAS_Cpp::Transp gpms_rhs_trans
    ,value_type beta
    ) const;

  /** \brief sym_lhs = alpha * op(mwo_rhs') * M * op(mwo_rhs).
    *
    * The default operation is based on <tt>this->Vp_StMtV(...)</tt> and assumes
    * that the matrix is symmetric.  Of course, a more efficient implementation
    * is often needed and the sublcass would like to override this.
    */
  virtual void Mp_StMtMtM(
    DMatrixSliceSym* sym_lhs, value_type alpha
    ,EMatRhsPlaceHolder dummy_place_holder
    ,const MatrixOpSerial& mwo_rhs, BLAS_Cpp::Transp mwo_rhs_trans
    ,value_type beta
    ) const;

  /** @name Overridden from MatrixSymOp */
  //@{

  /// Must be overridden to call <tt>MatrixOpSerial::space_rows()</tt>
  const VectorSpace& space_rows() const;

  /// symwo_lhs = alpha * op(gpms_rhs') * M * op(gpms_rhs) + beta * sym_lhs.
   void Mp_StPtMtP(
    MatrixSymOp* symwo_lhs, value_type alpha
    ,EMatRhsPlaceHolder dummy_place_holder
    ,const GenPermMatrixSlice& gpms_rhs, BLAS_Cpp::Transp gpms_rhs_trans
    ,value_type beta
    ) const;

  /// symwo_lhs = alpha * op(mwo_rhs') * M * op(mwo_rhs).
  void Mp_StMtMtM(
    MatrixSymOp* symwo_lhs, value_type alpha
    ,EMatRhsPlaceHolder dummy_place_holder
    ,const MatrixOp& mwo_rhs, BLAS_Cpp::Transp mwo_rhs_trans
    ,value_type beta
    ) const;

  //@}

};	// end class MatrixSymOpSerial

/** \defgroup MatrixSymWithOpSerial_funcs Inline nonmeber functions for <tt>MatrixSymOpSerial</tt>.
  */
//@{

inline
/// sym_lhs = alpha * op(gpms_rhs') * M * op(gpms_rhs) + beta * sym_lhs.
void Mp_StPtMtP(
  DMatrixSliceSym* sym_lhs, value_type alpha
  ,MatrixSymOpSerial::EMatRhsPlaceHolder dummy_place_holder
  ,const MatrixSymOpSerial& M
  ,const GenPermMatrixSlice& gpms_rhs, BLAS_Cpp::Transp gpms_rhs_trans
  ,value_type beta = 1.0
  )
{
  M.Mp_StPtMtP(sym_lhs,alpha,dummy_place_holder,gpms_rhs,gpms_rhs_trans,beta);
}

inline
/// sym_lhs = alpha * op(mwo_rhs') * M * op(mwo_rhs) + beta * sym_lhs
void Mp_StMtMtM(
  DMatrixSliceSym* sym_lhs, value_type alpha
  ,MatrixSymOpSerial::EMatRhsPlaceHolder dummy_place_holder
  ,const MatrixSymOpSerial& M
  ,const MatrixOpSerial& mwo_rhs, BLAS_Cpp::Transp mwo_rhs_trans
  ,value_type beta = 1.0
  )
{
  M.Mp_StMtMtM(sym_lhs,alpha,dummy_place_holder,mwo_rhs,mwo_rhs_trans,beta);
}

//@}

}	// end namespace AbstractLinAlgPack 

#endif	// SLAP_MATRIX_SYM_WITH_OP_SERIAL_H
