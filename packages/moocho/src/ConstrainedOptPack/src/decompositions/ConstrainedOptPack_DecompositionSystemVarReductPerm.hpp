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

#ifndef DECOMPOSITION_SYSTEM_VAR_REDUCT_PERM_H
#define DECOMPOSITION_SYSTEM_VAR_REDUCT_PERM_H

#include <stdexcept>

#include "ConstrainedOptPack_DecompositionSystemVarReduct.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"

namespace ConstrainedOptPack {

/** \brief Specialization interface of \c DecompositonSystem that allows basis permutations.
 *
 * ToDo: Finish documentation!
 */
class DecompositionSystemVarReductPerm : public DecompositionSystemVarReduct {
public:

  /** @name Public types */
  //@{

  /** \brief . */
  typedef Teuchos::RCP<
    const Teuchos::AbstractFactory<Permutation> >         perm_fcty_ptr_t;

  //@}

  /** @name Constructors / initializers */
  //@{

  /** \brief . */
  DecompositionSystemVarReductPerm(
    EExplicitImplicit     D_imp    = MAT_IMP_AUTO
    ,EExplicitImplicit    Uz_imp   = MAT_IMP_AUTO
    )
    :DecompositionSystemVarReduct(D_imp,Uz_imp)
  {}

  //@}

  /** @name Permutation factories */
  //@{

  /** \brief . */
  virtual const perm_fcty_ptr_t   factory_P_var() const = 0;
  /** \brief . */
  virtual const perm_fcty_ptr_t   factory_P_equ() const = 0;

  //@}

  /** @name Setting or selecting a new decomposition */
  //@{

  /** \brief Query to see if a current basis is already selected.
   */
  virtual bool has_basis() const = 0;

  /** \brief Client selects the basis for <tt>Gc(:,con_decomp)'</tt>.
   *
   * ToDo: Finish documentation!
   */
  virtual void set_decomp(
    std::ostream          *out
    ,EOutputLevel         olevel
    ,ERunTests            test_what
    ,const Permutation    &P_var
    ,const Range1D        &var_dep
    ,const Permutation    *P_equ
    ,const Range1D        *equ_decomp
    ,const MatrixOp       &Gc
    ,MatrixOp             *Z
    ,MatrixOp             *Y
    ,MatrixOpNonsing      *R
    ,MatrixOp             *Uz
    ,MatrixOp             *Uy
    ,EMatRelations        mat_rel = MATRICES_INDEP_IMPS
    ) = 0;
  
  /** \brief Client asks decompostion system object to select the basis for <tt>Gc(:,con_decomp)'</tt>.
   *
   * ToDo: Finish documentation!
   */
  virtual void select_decomp(
    std::ostream              *out
    ,EOutputLevel             olevel
    ,ERunTests                test_what
    ,const Vector             *nu
    ,MatrixOp                 *Gc
    ,Permutation              *P_var
    ,Range1D                  *var_dep
    ,Permutation              *P_equ
    ,Range1D                  *equ_decomp
    ,MatrixOp                 *Z
    ,MatrixOp                 *Y
    ,MatrixOpNonsing          *R
    ,MatrixOp                 *Uz
    ,MatrixOp                 *Uy
    ,EMatRelations            mat_rel = MATRICES_INDEP_IMPS
    ) = 0;

  //@}
  
};	// end class DecompositionSystemVarReductPerm

}	// end namespace ConstrainedOptPack

#endif // DECOMPOSITION_SYSTEM_VAR_REDUCT_PERM_H
