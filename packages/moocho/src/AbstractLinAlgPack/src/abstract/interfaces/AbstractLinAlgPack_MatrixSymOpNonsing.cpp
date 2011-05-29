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

#include "AbstractLinAlgPack_MatrixSymOpNonsing.hpp"

namespace AbstractLinAlgPack {

MatrixSymOpNonsing::mat_mswons_mut_ptr_t
MatrixSymOpNonsing::clone_mswons()
{
  return Teuchos::null;
}

MatrixSymOpNonsing::mat_mswons_ptr_t
MatrixSymOpNonsing::clone_mswons() const
{
  return Teuchos::null;
}

// Overridden from MatrixOp

MatrixSymOpNonsing::mat_mut_ptr_t
MatrixSymOpNonsing::clone()
{
  return clone_mswons();
}

MatrixSymOpNonsing::mat_ptr_t
MatrixSymOpNonsing::clone() const
{
  return clone_mswons();
}

// Overridden from MatrixNonsing

MatrixSymOpNonsing::mat_mns_mut_ptr_t
MatrixSymOpNonsing::clone_mns()
{
  return clone_mswons();
}

MatrixSymOpNonsing::mat_mns_ptr_t
MatrixSymOpNonsing::clone_mns() const
{
  return clone_mswons();
}

// Overridden from MatrixSymOp

MatrixSymOpNonsing::mat_mswo_mut_ptr_t
MatrixSymOpNonsing::clone_mswo()
{
  return clone_mswons();
}

MatrixSymOpNonsing::mat_mswo_ptr_t
MatrixSymOpNonsing::clone_mswo() const
{
  return clone_mswons();
}

// Overridden from MatrixSymNonsing

MatrixSymOpNonsing::mat_msns_mut_ptr_t
MatrixSymOpNonsing::clone_msns()
{
  return clone_mswons();
}

MatrixSymOpNonsing::mat_msns_ptr_t
MatrixSymOpNonsing::clone_msns() const
{
  return clone_mswons();
}

// Overridden from MatrixOpNonsing

MatrixSymOpNonsing::mat_mwons_mut_ptr_t
MatrixSymOpNonsing::clone_mwons()
{
  return clone_mswons();
}

MatrixSymOpNonsing::mat_mwons_ptr_t
MatrixSymOpNonsing::clone_mwons() const
{
  return clone_mswons();
}

}	// end namespace AbstractLinAlgPack
