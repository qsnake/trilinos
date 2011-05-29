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

#include "NLPInterfacePack_test_basis_system.hpp"
#include "NLPInterfacePack_NLPFirstOrder.hpp"
#include "AbstractLinAlgPack_BasisSystem.hpp"
#include "AbstractLinAlgPack_BasisSystemTester.hpp"
#include "AbstractLinAlgPack_BasisSystemTesterSetOptions.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsing.hpp"

bool NLPInterfacePack::test_basis_system(
   NLPFirstOrder                                 *nlp
  ,BasisSystem                                  *basis_sys
  ,OptionsFromStreamPack::OptionsFromStream     *options
  ,std::ostream                                 *out
  )
{
  namespace mmp = MemMngPack;

  const index_type
    n = nlp->n(),
    m = nlp->m();

  // Create the matrices Gc and Gh
  NLPFirstOrder::mat_fcty_ptr_t::element_type::obj_ptr_t
    Gc = ( m  ? nlp->factory_Gc()->create() : Teuchos::null );
  
  // Compute the matrices at xinit
  const Vector
    &xo = nlp->xinit();
  if(m)
    nlp->set_Gc(Gc.get());
  if(m)
    nlp->calc_Gc(xo);

  // Create the matrices C and D
  BasisSystem::mat_nonsing_fcty_ptr_t::element_type::obj_ptr_t
    C = ( m ? basis_sys->factory_C()->create() : Teuchos::null);
  BasisSystem::mat_fcty_ptr_t::element_type::obj_ptr_t
    D = ( m && n > m && basis_sys->factory_C().get() ? basis_sys->factory_C()->create() : Teuchos::null);
  BasisSystem::mat_fcty_ptr_t::element_type::obj_ptr_t
    GcUP = ( m && n > m && basis_sys->factory_GcUP().get()  ? basis_sys->factory_GcUP()->create() : Teuchos::null);

  // Initialize C and D with basis_sys
  basis_sys->update_basis(
    *Gc
    ,C.get()
    ,D.get()
    ,GcUP.get()
    );

  // Test the basis and basis system objects.
  BasisSystemTester
    basis_sys_tester;
  if(options) {
    BasisSystemTesterSetOptions
      opt_setter(&basis_sys_tester);
    opt_setter.set_options(*options);
  }
  const bool result = basis_sys_tester.test_basis_system(
    *basis_sys
    ,Gc.get()
    ,C.get()
    ,NULL    // Create the N matrix internally
    ,D.get()
    ,GcUP.get()
    ,out
    );

  return result;
}
