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

#ifndef ABSTRACT_LIN_ALG_PACK_BASIS_PERM_DIRECT_SPARSE_SYSTEM_H
#define ABSTRACT_LIN_ALG_PACK_BASIS_PERM_DIRECT_SPARSE_SYSTEM_H

#include "AbstractLinAlgPack_DirectSparseSolver.hpp"
#include "AbstractLinAlgPack_BasisSystemPerm.hpp"
#include "DenseLinAlgPack_IVector.hpp"

namespace AbstractLinAlgPack {

/** \brief Permutatble basis system subclass that uses a direct sparse solver.
 *
 * This current implementation only allows undecomposed general inequality
 * constraints.
 *
 * ToDo: Finish documentation!
 */
class BasisSystemPermDirectSparse
  : public AbstractLinAlgPack::BasisSystemPerm
{
public:

  /** @name Public types */
  //@{

  /** \brief . */
  typedef Teuchos::RCP<DirectSparseSolver>   direct_solver_ptr_t;

  //@}

  /** @name Constructors / initializers */
  //@{

  /// Calls <tt>this->initialize()</tt>
  BasisSystemPermDirectSparse(
    const direct_solver_ptr_t&   direct_solver = Teuchos::null
    );
  
  /// Initialize given a direct sparse solver object.
  void initialize(
    const direct_solver_ptr_t&   direct_solver
    );

  //@}

  /** @name Overridden from BasisSystem */
  //@{

  /** \brief . */
  const mat_nonsing_fcty_ptr_t factory_C() const;
  /** \brief . */
  const mat_fcty_ptr_t factory_D() const;
  /** \brief . */
  const mat_fcty_ptr_t factory_GcUP() const;
  /** \brief . */
  Range1D var_dep() const;
  /** \brief . */
  Range1D var_indep() const;
  /** \brief . */
  Range1D equ_decomp() const;
  /** \brief . */
  Range1D equ_undecomp() const;
  /** \brief . */
  void update_basis(
    const MatrixOp          &Gc
    ,MatrixOpNonsing        *C
    ,MatrixOp               *D
    ,MatrixOp               *GcUP
    ,EMatRelations          mat_rel
    ,std::ostream           *out
    ) const;

  //@}

  /** @name Overridded from BasisSystemPerm */
  //@{

  /** \brief . */
  const perm_fcty_ptr_t   factory_P_var() const;
  /** \brief . */
  const perm_fcty_ptr_t   factory_P_equ() const;
  /** \brief . */
  const perm_fcty_ptr_t   factory_P_inequ() const;
  /** \brief . */
  void set_basis(
    const Permutation          &P_var
    ,const Range1D             &var_dep
    ,const Permutation         *P_equ
    ,const Range1D             *equ_decomp
    ,const MatrixOp            &Gc
    ,MatrixOpNonsing           *C
    ,MatrixOp                  *D
    ,MatrixOp                  *GcUP
    ,EMatRelations             mat_rel
    ,std::ostream              *out
    );
  /** \brief . */
  void select_basis(
    const Vector               *nu
    ,MatrixOp                  *Gc
    ,Permutation               *P_var
    ,Range1D                   *var_dep
    ,Permutation               *P_equ
    ,Range1D                   *equ_decomp
    ,MatrixOpNonsing           *C
    ,MatrixOp                  *D
    ,MatrixOp                  *GcUP
    ,EMatRelations             mat_rel
    ,std::ostream              *out
    );
  
  //@}

private:

  // ///////////////////////////////
  // Private data members

  direct_solver_ptr_t   direct_solver_;
  size_type             n_;
  size_type	          m_;
  size_type             r_;
  size_type             Gc_nz_;
  Range1D               init_var_rng_;
  IVector               init_var_inv_perm_;  // If init_var_rng is full range then this is ignored
    Range1D               init_equ_rng_;
  IVector               init_equ_inv_perm_;  // If init_equ_rng is full range then this is ignored
  Range1D               var_dep_;       // used by factor()
  Range1D               var_indep_;     // used by factor()
    Range1D               equ_decomp_;    // used by factor()
    Range1D               equ_undecomp_;  // used by factor()

  // ///////////////////////////////
  // Private member functions

  /** \brief . */
  Teuchos::RCP<DirectSparseSolver::BasisMatrix>
  get_basis_matrix( MatrixOpNonsingAggr &C_aggr ) const;

  /** \brief . */
  void set_A_mctse(
    size_type                    n
    ,size_type                   m
    ,const MatrixPermAggr        &Gc_pa
    ,MatrixConvertToSparseEncap  *A_mctse
    ) const;

  /** \brief . */
  void update_basis_and_auxiliary_matrices(
    const MatrixOp& Gc
    ,const Teuchos::RCP<DirectSparseSolver::BasisMatrix>& C_bm
    ,MatrixOpNonsingAggr *C_aggr
    ,MatrixOp* D, MatrixOp* GcUP
    ) const;

  /** \brief . */
  void do_some_basis_stuff(
    const MatrixOp& Gc
    ,const Range1D& var_dep, const Range1D& equ_decomp
    ,const Teuchos::RCP<DirectSparseSolver::BasisMatrix>& C_bm
    ,MatrixOpNonsingAggr *C_aggr
    ,MatrixOp* D, MatrixOp* GcUP
    );


}; // end class BasisSystemPermDirectSparse

}  // end namespace AbstractLinAlgPack

#endif // ABSTRACT_LIN_ALG_PACK_BASIS_PERM_DIRECT_SPARSE_SYSTEM_H
