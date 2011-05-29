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

#ifndef	ALAP_DIRECT_SPARSE_SOLVER_IMP_H
#define ALAP_DIRECT_SPARSE_SOLVER_IMP_H

#include "AbstractLinAlgPack_DirectSparseSolver.hpp"
#include "AbstractLinAlgPack_VectorSpaceSerial.hpp"

namespace AbstractLinAlgPack {

class DirectSparseSolverImp;

/** \brief Implementation node class for \c DirectSparseSolver that takes
 * care of the memory management details.
 *
 * ToDo: Finish documentation!
 *
 * Subclasses must override the following methods which are pure
 * virtual and are not implemented here: \c basis_matrix_factory(),
 * \c estimated_fillin_ratio(), \c create_fact_struc(),
 * \c create_fact_nonzeros(), \c imp_analyze_and_factor() and
 * \c imp_factor().
 */
class DirectSparseSolverImp : public DirectSparseSolver {
public:

  /** @name Protected types */
  //@{

  /** \brief Abstract class for objects that represent the factorization
    * nonzeros of a particular matrix.
    *
    * This storage for the nonzeros can be reused over and over again
    * for factored matrices with the same structure.
    */
  class FactorizationNonzeros {
  public:
    /** \brief . */
    virtual ~FactorizationNonzeros() {}
  };

  /** \brief Implementation node subclass that combines factorization structure and
   * factorization nonzeros into a single basis matrix object.
   *
   * The only methods that subclasses must override for a complete
   * basis matrix object are:
   * 
   * <tt>AbstractLinAlgPack::MatrixNonsing::V_InvMtV()</tt>,
   * and  <tt>create_matrix()</tt>.
   *
   * Note that is if very important that subclasses do not maintain their
   * own copies of smart reference counting pointer objects to the
   * factorization structure or factorization nonzeros.  Instead, the 
   * methods \c get_fact_struc() and \c get_fact_nonzeros() should
   * be called inside of the implementation of the the \c V_InvMtV()
   * method and then \c dynamic_cast<> used to get at the concreate
   * objects.
   */
  class BasisMatrixImp : public BasisMatrix {
  public:

    /** @name Public types */
    //@{

    /** \brief . */
    typedef Teuchos::RCP<FactorizationNonzeros> fact_nonzeros_ptr_t;

    //@}

    /** @name Access */
    //@{

    /** \brief Return a reference to a smart pointer to the object that represents
     * the factorization nonzeros.
     *
     * Returning a reference to a \c RCP<> object verses returning
     * a \c RCP<> object itself is critical so that we can rely on
     * \c RCP<>::count() to tell us how many clients have a reference
     * to this object.
     */
    virtual const fact_nonzeros_ptr_t&  get_fact_nonzeros() const;

    //@}

    /** @name Overridden from MatrixBase */
    //@{

    /** \brief . */
    const VectorSpace& space_cols() const;
    /** \brief . */
    const VectorSpace& space_rows() const;
    /** \brief . */
    size_type rows() const;
    /** \brief . */
    size_type cols() const;

    //@}

    /** @name Overridden from MatrixNonsinguar */
    //@{

    /** \brief . */
    mat_mns_mut_ptr_t clone_mns();

    //@}

    /** @name Overridden from BasisMatrix */
    //@{

    /** \brief . */
    virtual const fact_struc_ptr_t&  get_fact_struc() const;

    //@}

  protected:

    /** @name Constructors/initailizers */
    //@{

    /** \brief Default initializers to uninitialized.
     */
    BasisMatrixImp();

    /** \brief Calls <tt>this->initialize()</tt>
     */
    BasisMatrixImp(
      size_type                      dim
      ,const fact_struc_ptr_t        &fact_struc
      ,const fact_nonzeros_ptr_t     &fact_nonzeros
      );

    /** \brief Initialize given initialized factorization structure and factorization nonzeros objects.
     */
    virtual void initialize(
      size_type                      dim
      ,const fact_struc_ptr_t        &fact_struc
      ,const fact_nonzeros_ptr_t     &fact_nonzeros
      );

    /** \brief Make uninitialized.
     *
     * Postconditions:<ul>
     * <li> <tt>this->dim() == 0</tt>
     * <li> <tt>this->get_fact_struc().get() == NULL</tt>
     * <li> <tt>this->get_fact_nonzeros().get() == NULL</tt>
     * </ul>
     */
    void set_uninitialized();

    //@}

    /** @name Factory methods to be overridden by subclasses */
    //@{

    /// Called by \c this->clone-mns(). 
    virtual Teuchos::RCP<BasisMatrixImp> create_matrix() const = 0;

    //@}

  private:

    /** \brief Allow only DirectSparseSolverImp to initialize objects of this type.
     * Important !!!!! Even though DirectSparseSolverImp has access to these
     * private members it is strictly not to access them directly !!!!
     */
    friend class DirectSparseSolverImp;

#ifdef DOXYGEN_COMPILE
    FactorizationStructure    *fact_struc;
    FactorizationNonzeros     *fact_nonzeros;
#else
    size_type                 dim_;
    fact_struc_ptr_t          fact_struc_;
    fact_nonzeros_ptr_t       fact_nonzeros_;
    VectorSpaceSerial         vec_space_;
#endif

  }; // end class BasisMatrixImp

  //@}

  /** @name Overridden from DirectSparseSolver */
  //@{

  /** \brief . */
  void analyze_and_factor(
    const AbstractLinAlgPack::MatrixConvertToSparse   &A
    ,DenseLinAlgPack::IVector                            *row_perm
    ,DenseLinAlgPack::IVector                            *col_perm
    ,size_type                                      *rank
    ,BasisMatrix                                    *basis_matrix
    ,std::ostream                                   *out
    );
  /** \brief . */
  void factor(
    const AbstractLinAlgPack::MatrixConvertToSparse   &A
    ,BasisMatrix                                    *basis_matrix
    ,const BasisMatrix::fact_struc_ptr_t            &fact_struc
    ,std::ostream                                   *out
    );
  /** \brief . */
  const BasisMatrix::fact_struc_ptr_t& get_fact_struc() const;
  /** \brief . */
  void set_uninitialized();

  //@}

protected:

  /** @name Protected pure virtual methods to be overridden by concrete direct solver subclasses */
  //@{

  /** \brief Create a new, uninitialized \c FactorizationStructure object.
   */
  virtual const Teuchos::RCP<FactorizationStructure> create_fact_struc() const = 0;

  /** \brief Create a new, uninitialized \c FactorizationNonzeros object.
   */
  virtual const Teuchos::RCP<FactorizationNonzeros> create_fact_nonzeros() const = 0;

  /** \brief Called to implement the \c analyze_and_factor() without having to worry about
   * memory mangagment details.
   *
   * ToDo: Finish documentation!
   */	
  virtual void imp_analyze_and_factor(
    const AbstractLinAlgPack::MatrixConvertToSparse   &A
    ,FactorizationStructure                         *fact_struc
    ,FactorizationNonzeros                          *fact_nonzeros
    ,DenseLinAlgPack::IVector                            *row_perm
    ,DenseLinAlgPack::IVector                            *col_perm
    ,size_type                                      *rank
    ,std::ostream                                   *out            = NULL
    ) = 0;

  /** \brief Called to implement the \c analyze_and_factor() without having to worry about
   * memory mangagment details.
   *
   * ToDo: Finish documentation!
   */	
  virtual void imp_factor(
    const AbstractLinAlgPack::MatrixConvertToSparse   &A
    ,const FactorizationStructure                   &fact_struc
    ,FactorizationNonzeros                          *fact_nonzeros
    ,std::ostream                                   *out            = NULL
    ) = 0;

  //@}

private:

#ifdef DOXYGEN_COMPILE
    FactorizationStructure           *fact_struc;
#else
    BasisMatrix::fact_struc_ptr_t    fact_struc_;
    size_type                        rank_;
#endif

};	// end class DirectSparseSolverImp 

}	// end namespace AbstractLinAlgPack 

#endif	// ALAP_DIRECT_SPARSE_SOLVER_IMP_H
