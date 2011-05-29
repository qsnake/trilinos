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

#include "Moocho_ConfigDefs.hpp"


#ifdef HAVE_MOOCHO_MA28


#include "AbstractLinAlgPack_MA28Solver.hpp"

// Initialize static variables
MA28_Cpp::MA28Solver* MA28_Cpp::MA28Solver::curr_solver_ = 0;

// Initialize the references to the ma28 common blocks defined in the fortran code.
MA28_Cpp::MA28CommonBlockReferences MA28_Cpp::MA28Solver::ma28_common_blocks_(
    MA28_CppDecl::ma28ed_cb,
    MA28_CppDecl::ma28fd_cb,
    MA28_CppDecl::ma28gd_cb,
    MA28_CppDecl::ma28hd_cb,
    MA28_CppDecl::ma30ed_cb,
    MA28_CppDecl::ma30fd_cb,
    MA28_CppDecl::ma30gd_cb,
    MA28_CppDecl::ma30hd_cb,
    MA28_CppDecl::ma30id_cb,
    MA28_CppDecl::mc23bd_cb
  );

// Save the default values of the ma28 common block variables.
// By setting the references to the common block variables first we are guarented that
// references will be bound to the proper memory locations.  The only problem with this is that
// how can I be sure that the initializations for the fortran BLOCK DATA units have
// be performed before these global initializations are carried out.  I need to look into this in
// different platforms to ensure that this will work.  Otherwise I will just have to set the
// default values myself.
MA28_Cpp::MA28CommonBlockStorage MA28_Cpp::MA28Solver::default_common_blocks_(
  MA28Solver::ma28_common_blocks_);

// ///////////////////////////////
// Member functions

MA28_Cpp::MA28Solver::MA28Solver()
  : common_blocks_(default_common_blocks_), changed_(false)
{
  mp(0);	// We can't print to standard fortran streams from C++
  lp(0);
}

MA28_Cpp::MA28Solver::MA28Solver(const MA28Solver& s)
  : common_blocks_(s.common_blocks_), changed_(false)
{}

void MA28_Cpp::MA28Solver::set_common_block_data() {
  // Copy to the ma28 ma28 common blocks if this is not the current solver.
  if( (this != curr_solver_ )|| ( (this == curr_solver_) && changed_ ) )
    ma28_common_blocks_ = common_blocks_;
  // Make this solver the current solver
  curr_solver_ = this;
  changed_ = false;
}

void MA28_Cpp::MA28Solver::get_common_block_data() {
  // You must save the common block data back in order to get the return parameters.
  // Later I should differentiate what is control and what is return but for
  // now this is the easiest thing to do.
  common_blocks_ = ma28_common_blocks_;
}


#endif // HAVE_MOOCHO_MA28
