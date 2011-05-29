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

#ifndef ABSTRACT_LIN_ALG_PACK_BASIS_PERM_SYSTEM_H
#define ABSTRACT_LIN_ALG_PACK_BASIS_PERM_SYSTEM_H

#include "AbstractLinAlgPack_BasisSystem.hpp"

namespace AbstractLinAlgPack {

/** \brief Interface for setting and selecting a basis from the Jacobian
 * from a set of equations.
 *
 * ToDo: Finish documentation!
 */
class BasisSystemPerm : public BasisSystem {
public:

  /** @name Public types */
  //@{

  /** \brief . */
  typedef Teuchos::RCP<
    const Teuchos::AbstractFactory<Permutation> >   perm_fcty_ptr_t;

  //@}


  /** \brief Required constructor (calls <tt>initialize()</tt>).
   */
  BasisSystemPerm(
    const mat_sym_fcty_ptr_t             &factory_transDtD
    ,const mat_sym_nonsing_fcty_ptr_t    &factory_S
    )
    :BasisSystem(factory_transDtD,factory_S)
      {}

  /** @name Permutation factories */
  //@{

  /** \brief . */
  virtual const perm_fcty_ptr_t   factory_P_var() const = 0;
  /** \brief . */
  virtual const perm_fcty_ptr_t   factory_P_equ() const = 0;

  //@}

  /** @name Basis selection / manipulation */
  //@{

  /** \brief Factor a basis selected by the client.
   *
   * ToDo: Finish documentation!
   */
  virtual void set_basis(
    const Permutation          &P_var
    ,const Range1D             &var_dep
    ,const Permutation         *P_equ
    ,const Range1D             *equ_decomp
    ,const MatrixOp            &Gc
    ,MatrixOpNonsing           *C
    ,MatrixOp                  *D
    ,MatrixOp                  *GcUP
    ,EMatRelations             mat_rel = MATRICES_INDEP_IMPS
    ,std::ostream              *out    = NULL
    ) = 0;

  /** \brief Select a basis.
   *
   * ToDo: Finish documentation!
   */
  virtual void select_basis(
    const Vector               *nu
    ,MatrixOp                  *Gc
    ,Permutation               *P_var
    ,Range1D                   *var_dep
    ,Permutation               *P_equ
    ,Range1D                   *equ_decomp
    ,MatrixOpNonsing           *C
    ,MatrixOp                  *D
    ,MatrixOp                  *GcUP
    ,EMatRelations             mat_rel = MATRICES_INDEP_IMPS
    ,std::ostream              *out    = NULL
    ) = 0;
  
  //@}

private:
  // not defined and not to be called
  BasisSystemPerm();

}; // end class BasisSystemPerm

}  // end namespace AbstractLinAlgPack

#endif // ABSTRACT_LIN_ALG_PACK_BASIS_PERM_SYSTEM_H
