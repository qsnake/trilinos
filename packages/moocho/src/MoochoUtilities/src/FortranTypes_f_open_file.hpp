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

#ifndef F_OPEN_FILE_H
#define F_OPEN_FILE_H

#include "Moocho_ConfigDefs.hpp"
#include "Teuchos_F77_wrappers.h"

namespace FortranTypes {

/** @name Open a Fortran file.
  */
//@{

/** \brief . */
enum EOpenStatus { OPEN_OLD = 0, OPEN_NEW = 1, OPEN_SCRATCH = 2
  , OPEN_UNKNOWN = 3 };
/** \brief . */
enum EOpenForm { OPEN_FORMATTED = 0, OPEN_UNFORMATTED = 1 };
/** \brief . */
enum EOpenBlank { OPEN_NULL = 0, OPEN_ZERO = 1 };
/** \brief . */
enum EOpenAccess { OPEN_SEQUENTIAL = 0, OPEN_DIRECT = 1 };

/** Open a Fortran file given its name and unit number.
  *
  * If successful #iunit# is returned for the opened file.
  *
  * The standard options to Fortran OPEN(...) are included.  The
  * only mandatory option to set is for the file name in FILE.
  *
  * The length of the file name must be <= 100 characters long.
  *
  * Note that all of the options are optional but if
  * you set access == OPEN_DIRECT you must set recl to some
  * value greater than zero.  See Metcalf, 1990, p. 127.
  *
  * If #file# is not a valid ASCII string then the exception
  * InvalidFileNameException will be thrown
  *
  * If the file could not be opened for some reason then
  * the exception OpenException will be thrown.
  */
void f_open_file( const f_int iunit, const char file[]
  , EOpenStatus status = OPEN_UNKNOWN, EOpenForm form = OPEN_FORMATTED
  , EOpenBlank blank = OPEN_NULL, EOpenAccess access = OPEN_SEQUENTIAL
  , f_int recl = -1 );

/** Close a Fortran file given its unit number.
  *
  * If successful #iunit# is returned for the opened file.
  *
  * The standard options to Fortran OPEN(...) are included.  The
  * only mandatory option to set is for the file name in FILE.
  *
  * The length of the file name must be <= 100 characters long.
  *
  * Note that all of the options are optional but if
  * you set access == OPEN_DIRECT you must set recl to some
  * value greater than zero.  See Metcalf, 1990, p. 127.
  *
  * If #file# is not a valid ASCII string then the exception
  * InvalidFileNameException will be thrown
  *
  * If the file could not be opened for some reason then
  * the exception OpenException will be thrown.
  */
void f_close_file( const f_int iunit, bool keep = true );

//@}

/// Thrown if the file name is not a valid ASCII string
class InvalidFileNameException : public std::logic_error
{public: InvalidFileNameException(const std::string& what_arg) : std::logic_error(what_arg) {}};

/// Thrown if the open operation fails
class OpenException : public std::logic_error
{public: OpenException(const std::string& what_arg) : std::logic_error(what_arg) {}};

}	// end namespace FortranTypes 

#endif // F_OPEN_FILE_H
