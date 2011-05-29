
/*! @HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */


/*! @file CTrilinos_exceptions.hpp
 * @brief Defines exceptions thrown by CTrilinos. */


#ifndef CTRILINOS_EXCEPTIONS_HPP
#define CTRILINOS_EXCEPTIONS_HPP


#include "CTrilinos_config.h"


#include "Teuchos_Exceptions.hpp"


namespace CTrilinos {


/*! exception indicating wrong object type encountered */
class CTrilinosTypeMismatchError : public Teuchos::ExceptionBase
{
  public:
    CTrilinosTypeMismatchError(const std::string& what_arg) : Teuchos::ExceptionBase(what_arg) {}
};

/*! exception indicating wrong object type encountered */
class CTrilinosInvalidTypeError : public Teuchos::ExceptionBase
{
  public:
    CTrilinosInvalidTypeError(const std::string& what_arg) : Teuchos::ExceptionBase(what_arg) {}
};

/*! exception indicating invalid cast attempted */
class CTrilinosConstCastError : public Teuchos::ExceptionBase
{
  public:
    CTrilinosConstCastError(const std::string& what_arg) : Teuchos::ExceptionBase(what_arg) {}
};

/*! exception indicating wrong object table accessed */
class CTrilinosWrongTableError : public Teuchos::ExceptionBase
{
  public:
    CTrilinosWrongTableError(const std::string& what_arg) : Teuchos::ExceptionBase(what_arg) {}
};

/*! exception indicating some other error */
class CTrilinosMiscException : public Teuchos::ExceptionBase
{
  public:
    CTrilinosMiscException(const std::string& what_arg) : Teuchos::ExceptionBase(what_arg) {}
};


} // namespace CTrilinos


#endif

