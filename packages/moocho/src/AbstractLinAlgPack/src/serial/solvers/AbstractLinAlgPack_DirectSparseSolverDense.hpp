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

#ifndef	DIRECT_SPARSE_SOLVER_DENSE_H
#define DIRECT_SPARSE_SOLVER_DENSE_H

#include <valarray>
#include <vector>
#include <string>

#include "AbstractLinAlgPack_DirectSparseSolverImp.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "DenseLinAlgPack_IVector.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace AbstractLinAlgPack {

/** \brief Concreate sparse solver subclass that uses the dense LAPACK routines.
 *
 * ToDo: Finish documentation!
 */
class DirectSparseSolverDense : public DirectSparseSolverImp {
public:

  /** @name Constructors/initializers */
  //@{

  /** \brief Default constructor */
  DirectSparseSolverDense();

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

  /** \brief Implements the BasisMatrix object for Dense.
   */
  class BasisMatrixDense : public BasisMatrixImp {
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

  }; // end class BasisMatrixDense

  /** \brief Stores the factorization structure for Dense
   */
  class FactorizationStructureDense : public FactorizationStructure {
  public:
    friend class DirectSparseSolverDense;
    friend class BasisMatrixDense;
  private:
    FortranTypes::f_int      m_;         // Number of rows in A
    FortranTypes::f_int      n_;         // Number of columns in A
    FortranTypes::f_int      nz_;        // Number of nonzeros in A
    FortranTypes::f_int      rank_;      // Rank of the basis
    IVector                  col_perm_;  // First rank entries selects the basis of A
    IVector                  inv_col_perm_; // Inverse of col_perm_
    FactorizationStructureDense();
  }; // end class FactorizationStructureDense

  /** \brief Stores the factorization nonzeros for Dense
   */
  class FactorizationNonzerosDense : public FactorizationNonzeros {
  public:
    typedef FortranTypes::f_int f_int;
    friend class DirectSparseSolverDense;
    friend class BasisMatrixDense;
  private:
    DMatrix                          LU_;
    bool                               rect_analyze_and_factor_; // true for n > m analyze_and_factor()
    std::valarray<f_int>               ipiv_; // The permutation sent to xGETRS (identity if rect_analyze_and_factor_==true)
    IVector                            basis_perm_; // Only used if rect_analyze_and_factor_==true
  }; // end class FactorizationNonzerosDense

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

};	// end class DirectSparseSolverDense 

}	// end namespace AbstractLinAlgPack 

#endif	// DIRECT_SPARSE_SOLVER_DENSE_H
