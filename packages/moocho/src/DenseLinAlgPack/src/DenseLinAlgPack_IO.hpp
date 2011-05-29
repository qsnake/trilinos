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
//
// Includes all the headers for the DenseLinAlgPack input/output functions and
// stream operators for DVector, DVectorSlice, DMatrix, and DMatrixSlice.
// The utility function for "eating" comment lines is also included.

#ifndef LINALGPACK_IO_H
#define LINALGPACK_IO_H

#include "InputStreamHelperPack_EatInputComment.hpp"
#include "DenseLinAlgPack_DVectorIn.hpp"
#include "DenseLinAlgPack_DVectorOut.hpp"
#include "DenseLinAlgPack_DMatrixIn.hpp"
#include "DenseLinAlgPack_DMatrixOut.hpp"
#include "DenseLinAlgPack_InFormat.hpp"
#include "DenseLinAlgPack_OutFormat.hpp"

// Include namelookups for templated operator functions.  MS VS++ 5.0 standard
// nonconformance problem.

#include "DenseLinAlgPack_IO_NameLookups.hpp"

#endif // LINALGPACK_IO_H
