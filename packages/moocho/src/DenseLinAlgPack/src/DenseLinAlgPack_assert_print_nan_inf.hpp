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

#ifndef ASSERT_PRINT_NAN_INF_H
#define ASSERT_PRINT_NAN_INF_H

#include <stdexcept>

#include "DenseLinAlgPack_Types.hpp"

namespace DenseLinAlgPack {

/** \brief . */
class NaNInfException : public std::runtime_error
{public: NaNInfException(const std::string& what_arg) : std::runtime_error(what_arg) {}};

/** \brief . */
/* * This function asserts if a value_type scalare is a NaN or Inf and optionally
  * prints out these entires.
  * 
  * @param	val				[I]	Value the check
  * @param	name 			[I]	Name of the scale variable for output purposes
  * @param	throw_excpt		[I]	If true and is found to be a NaN or Inf
  * 							then a NaNInfException excetion is thrown after
  * 							any output.
  *	@param	out				[I/O]	If out==NULL then not output is produced.
  *									If out!=NULL and val is not
  *									NaN or Inf, then no output is produced.
  *									If out!=NULL and val is
  *									NaN or Inf then this will be printed before any
  *									execption is thrown.
  *									
  *	@return Returns true if val is not NaN or Inf.  If val
  *		is NaN or Inf then false will be returned unless a
  *		excetion NaNInfException was thrown (throw_except==true).
  */
bool assert_print_nan_inf( const value_type& val, char name[]
  , bool throw_excpt, std::ostream* out );

/** \brief . */
/* * This function asserts if a vector has any NaN or inf entries and optionally
  * prints out these entires.
  * 
  * @param	v 				[I]	DVector slice to check
  * @param	name 			[I]	Name of the vector for output purposes
  * @param	throw_excpt		[I]	If true and an entry is found to be a NaN or Inf
  * 							then a NaNInfException excetion is thrown after
  * 							any output.
  *	@param	out				[I/O]	If out==NULL then not output is produced.
  *									If out!=NULL and none of the entries is
  *									NaN or Inf, then no output is produced.
  *									If out!=NULL then any entries that are
  *									NaN or Inf will be printed before any
  *									execption is thrown.
  *									
  *	@return Returns true none of the entries are NaN or Inf.  If one of the
  *		entries is NaN or Inf then false will be returned unless an
  *		excetion was thrown (throw_except==true).
  */
bool assert_print_nan_inf( const DVectorSlice& v, char name[]
  , bool throw_excpt, std::ostream* out );

/** \brief . */
/* * This function asserts if a matrix has any NaN or inf entries and optionally
  * prints out these entires.
  * 
  * @param	m 				[I]	Matrix slice to check
  * @param	name 			[I]	Name of the matrix for output purposes
  * @param	throw_excpt		[I]	If true and an entry is found to be a NaN or Inf
  * 							then a NaNInfException excetion is thrown after
  * 							any output.
  *	@param	out				[I/O]	If out==NULL then not output is produced.
  *									If out!=NULL and none of the entries is
  *									NaN or Inf, then no output is produced.
  *									If out!=NULL then any entries that are
  *									NaN or Inf will be printed before any
  *									execption is thrown.
  *									
  *	@return Returns true none of the entries are NaN or Inf.  If one of the
  *		entries is NaN or Inf then false will be returned unless an
  *		excetion was thrown (throw_except==true).
  */
bool assert_print_nan_inf( const DMatrixSlice& m, char name[]
  , bool throw_excpt, std::ostream* out );

}	// end namespace DenseLinAlgPack

#endif // ASSERT_PRINT_NAN_INF_H
