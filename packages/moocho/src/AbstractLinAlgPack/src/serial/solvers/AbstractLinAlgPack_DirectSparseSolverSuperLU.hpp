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

#ifdef SPARSE_SOLVER_PACK_USE_SUPERLU

#ifndef	DIRECT_SPARSE_SOLVER_SUPERLU_H
#define DIRECT_SPARSE_SOLVER_SUPERLU_H

#include <valarray>
#include <vector>
#include <string>

#include "AbstractLinAlgPack_DirectSparseSolverImp.hpp"
#include "AbstractLinAlgPack_SuperLUSolver.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_IVector.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace AbstractLinAlgPack {

/** \brief Concreate sparse solver subclass that uses SuperLU.
 *
 * ToDo: Finish documentation!
 */
class DirectSparseSolverSuperLU : public DirectSparseSolverImp {
public:

  /** @name Control parameters */
  //@{

  // ToDo: Fill these in!

  //@}

  /** @name Constructors/initializers */
  //@{

  /** \brief Default constructor */
  DirectSparseSolverSuperLU();

  //@}

  /** @name Overridden from DirectSparseSolver */
  //@{

  /** \brief . */
  const basis_matrix_factory_ptr_t basis_matrix_factory() const;
  /** \brief . */
  void estimated_fillin_ratio( value_type estimated_fillin_ratio );

  //@}

protected:

  /** @name Protected types */
  //@{

  /** \brief Implements the BasisMatrix object for SuperLU.
   */
  class BasisMatrixSuperLU : public BasisMatrixImp {
  public:

    /** @name Overridden from BasisMatrixImp */
    //@{

    /** \brief . */
    Teuchos::RCP<BasisMatrixImp> create_matrix() const;
    /** \brief . */
    void V_InvMtV(
      VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
      ,const Vector& v_rhs2) const ;
    
    //@}

  }; // end class BasisMatrixSuperLU

  /** \brief Stores the factorization structure for SuperLU
   */
  class FactorizationStructureSuperLU : public FactorizationStructure {
  public:
    friend class DirectSparseSolverSuperLU;
    friend class BasisMatrixSuperLU;
  private:
    Teuchos::RCP<SuperLUPack::SuperLUSolver>
      superlu_solver_;
    Teuchos::RCP<SuperLUPack::SuperLUSolver::FactorizationStructure>
      fact_struct_;
    FactorizationStructureSuperLU();
  }; // end class FactorizationStructureSuperLU

  /** \brief Stores the factorization nonzeros for SuperLU
   */
  class FactorizationNonzerosSuperLU : public FactorizationNonzeros {
  public:
    friend class DirectSparseSolverSuperLU;
    friend class BasisMatrixSuperLU;
  private:
    Teuchos::RCP<SuperLUPack::SuperLUSolver::FactorizationNonzeros>
      fact_nonzeros_;
  }; // end class FactorizationNonzerosSuperLU

  //@}

  /** @name Overridden from DirectSparseSolverImp */
  //@{

  /** \brief . */
  const Teuchos::RCP<FactorizationStructure> create_fact_struc() const;
  /** \brief . */
  const Teuchos::RCP<FactorizationNonzeros> create_fact_nonzeros() const;
  /** \brief . */
  void imp_analyze_and_factor(
    const AbstractLinAlgPack::MatrixConvertToSparse   &A
    ,FactorizationStructure                         *fact_struc
    ,FactorizationNonzeros                          *fact_nonzeros
    ,DenseLinAlgPack::IVector                            *row_perm
    ,DenseLinAlgPack::IVector                            *col_perm
    ,size_type                                      *rank
    ,std::ostream                                   *out
    );
  /** \brief . */
  void imp_factor(
    const AbstractLinAlgPack::MatrixConvertToSparse   &A
    ,const FactorizationStructure                   &fact_struc
    ,FactorizationNonzeros                          *fact_nonzeros
    ,std::ostream                                   *out
    );

  //@}

private:

  // /////////////////////////////////
  // Private data members

  // ////////////////////////////////
  // Private member functions

};	// end class DirectSparseSolverSuperLU 

// ////////////////////////////////////////
// Inline members

}	// end namespace AbstractLinAlgPack 

#endif	// DIRECT_SPARSE_SOLVER_SUPERLU_H

#endif // SPARSE_SOLVER_PACK_USE_SUPERLU
