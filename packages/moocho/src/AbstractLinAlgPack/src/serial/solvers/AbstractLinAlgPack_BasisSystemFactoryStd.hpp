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

#ifndef SPARSE_SOLVER_PACK_BASIS_SYSTEM_FACTORY_H
#define SPARSE_SOLVER_PACK_BASIS_SYSTEM_FACTORY_H

#include "AbstractLinAlgPack_Types.hpp"
#include "AbstractLinAlgPack_BasisSystemFactory.hpp"
#include "Teuchos_AbstractFactory.hpp"

namespace AbstractLinAlgPack {

/** \brief Default implementation for <tt>BasisSystemPermDirectSparse</tt> obejcts
 * using <tt>DirectSparseSolver</tt> object.
 *
 * Several direct sparse solvers are supported by default.  These include:
 * <ul>
 * <li> DENSE (using LAPACK xGETRF())
 * <li> MA28
 * <li> MA48 (using MA28 for BasisSystemPerm::select_basis()) (not yet)
 * <li> SuperLU
 * </ul>
 *
 * These solvers are supported only if the proper macros are defined.
 * 
 * ToDo: Create a DirectSparseSolverFactory interface and use this
 * to allow clients to add new DirectSparseSolvers ...
 *
 */
class BasisSystemFactoryStd
  : public AbstractLinAlgPack::BasisSystemFactory
{
public:

  /** \brief . */
  BasisSystemFactoryStd(); // ToDo: Add arguments!

  /** @name Overridden from BasisSystemFactory */
  //@{

  /** \brief . */
  void set_options( const options_ptr_t& options );
  /** \brief . */
  const options_ptr_t& get_options() const;

  //@}

  /** @name Overridden from AbstractFactory */
  //@{

  /** \brief . */
  obj_ptr_t create() const;

  //@}

private:

  // ////////////////////////
  // Private types

  enum EDirectLinearSolverType { LA_DENSE, LA_MA28, LA_MA48, LA_SUPERLU };

  // ////////////////////////
  // Private data members

  mutable EDirectLinearSolverType  direct_linear_solver_type_;
  options_ptr_t                    options_;

  // ////////////////////////
  // Private member functions
  
  void read_options() const;

}; // end class BasisSystemFactoryStd

}  // end namespace AbstractLinAlgPack

#endif // SPARSE_SOLVER_PACK_BASIS_SYSTEM_FACTORY_H
