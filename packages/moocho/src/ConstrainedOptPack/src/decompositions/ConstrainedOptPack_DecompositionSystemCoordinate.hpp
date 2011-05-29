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

#ifndef DECOMPOSITION_SYSTEM_COORDINATE_H
#define DECOMPOSITION_SYSTEM_COORDINATE_H

#include "ConstrainedOptPack_DecompositionSystemVarReductImp.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace ConstrainedOptPack {

/** \brief Coordinate variable reduction subclass.
 *
 * This is the interface for the coordinate variable reduction decomposition
 * where:
 \verbatim

  Y = [ I ]   (class MatrixIdentConcatStd with MatrixZero)
      [ 0 ]

  R = Gc(:,con_decomp)'*Y = [ C   N ] * [ I ] = C
                                        [ 0 ]

  Uy = Gc(:,con_undecomp)'*Y = [ E  F ] * [ I ] = E
                                          [ 0 ]

 \endverbatim
 * The solution of the
 *
 * For now the copy constructor and the assignment operator are not defined.
 */
class DecompositionSystemCoordinate : public DecompositionSystemVarReductImp {
public:

  /** @name Constructors / initializers */
  //@{

  /** \brief . */
  DecompositionSystemCoordinate(
    const VectorSpace::space_ptr_t     &space_x          = Teuchos::null
    ,const VectorSpace::space_ptr_t    &space_c          = Teuchos::null
    ,const basis_sys_ptr_t             &basis_sys        = Teuchos::null
    ,const basis_sys_tester_ptr_t      &basis_sys_tester = Teuchos::null
    ,EExplicitImplicit                 D_imp             = MAT_IMP_AUTO
    ,EExplicitImplicit                 Uz_imp            = MAT_IMP_AUTO
    );

  //@}

  /** @name Overridden from DecompositionSystem */
  //@{

  /** \brief . */
  const mat_fcty_ptr_t factory_Y() const;
  /** \brief . */
  const mat_nonsing_fcty_ptr_t factory_R() const;
  /** \brief . */
  const mat_fcty_ptr_t factory_Uy() const;

  //@}

protected:

  /** @name Overridden from DecompositionSystemVarReductImp */
  //@{

  /** \brief . */
  mat_nonsing_fcty_ptr_t::element_type::obj_ptr_t	uninitialize_matrices(
    std::ostream                                       *out
    ,EOutputLevel                                      olevel
    ,MatrixOp                                          *Y
    ,MatrixOpNonsing                                   *R
    ,MatrixOp                                          *Uy
    ) const;
  /** \brief . */
  void initialize_matrices(
    std::ostream                                           *out
    ,EOutputLevel                                          olevel
    ,const mat_nonsing_fcty_ptr_t::element_type::obj_ptr_t &C
    ,const mat_fcty_ptr_t::element_type::obj_ptr_t         &D
    ,MatrixOp                                              *Y
    ,MatrixOpNonsing                                       *R
    ,MatrixOp                                              *Uy
    ,EMatRelations                                         mat_rel
    ) const;
  /** \brief . */
  void print_update_matrices(
    std::ostream& out, const std::string& leading_str ) const;

  //@}

};	// end class DecompositionSystemCoordinate

}	// end namespace ConstrainedOptPack

#endif	// DECOMPOSITION_SYSTEM_COORDINATE_H
