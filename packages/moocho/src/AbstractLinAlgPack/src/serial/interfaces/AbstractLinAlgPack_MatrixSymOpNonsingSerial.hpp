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

#ifndef SLAP_MATRIX_SYM_WITH_OP_NONSINGULAR_SERIAL_H
#define SLAP_MATRIX_SYM_WITH_OP_NONSINGULAR_SERIAL_H

#include "AbstractLinAlgPack_MatrixSymOpSerial.hpp"
#include "AbstractLinAlgPack_MatrixSymNonsingSerial.hpp"
#include "AbstractLinAlgPack_MatrixSymOpNonsing.hpp"

namespace AbstractLinAlgPack {

/** \brief Abstract base class for all serial polymorphic symmetric nonsingular matrices that
  * can be used to compute matrix-vector products and solve for linear systems relatively
  * efficiently.
  */
class MatrixSymOpNonsingSerial 
  : virtual public MatrixSymOpSerial
  , virtual public MatrixSymNonsingSerial
  , virtual public AbstractLinAlgPack::MatrixOpNonsing      // doxygen needs full name
  , virtual public AbstractLinAlgPack::MatrixSymOpNonsing   // ""
{
  // These overrides are needed for MS Visual C++ 2008
  virtual mat_mut_ptr_t clone() { return MatrixOp::clone(); }
  virtual mat_ptr_t clone() const { return MatrixOp::clone(); }
  virtual mat_mns_mut_ptr_t clone_mns() { return MatrixNonsing::clone_mns(); }
  virtual mat_mns_ptr_t clone_mns() const { return MatrixNonsing::clone_mns(); }
};

}	// end namespace AbstractLinAlgPack

#endif	// SLAP_MATRIX_SYM_WITH_OP_NONSINGULAR_SERIAL_H
