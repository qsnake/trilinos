/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
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
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#ifndef SUNDANCE_EXCEPTIONS_H
#define SUNDANCE_EXCEPTIONS_H

#include "SundanceDefs.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_TestForException.hpp"
#include <stdexcept>

#ifndef DOXYGEN_DEVELOPER_ONLY

//bvbw for backard compatibility reasons
//     I could not get this to work with ifdefs hence the hack

#define SUNDANCE_ERROR7(msg) \
{ \
  TestForException_break(); \
  TeuchosOStringStream omsg; \
	omsg << __FILE__ << ":" << __LINE__ << ": " \
       << ": " << msg; \
  throw Sundance::RuntimeError(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg)); \
}

#define SUNDANCE_ERROR(msg) \
{ \
  TeuchosOStringStream omsg; \
	omsg << __FILE__ << ":" << __LINE__ << ": " \
       << ": " << msg; \
  const std::string &omsgstr = omsg.str(); \
  TestForException_break(omsgstr); \
  throw Sundance::RuntimeError(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg)); \
}


#define SUNDANCE_TRACE(e) \
{ \
  TeuchosOStringStream omsg; \
	omsg << e.what() << std::endl \
  << "caught in " << __FILE__ << ":" << __LINE__ << std::endl ; \
  throw Sundance::RuntimeError(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg)); \
}

#define SUNDANCE_TRACE_MSG(e, msg)                      \
{ \
  TeuchosOStringStream omsg; \
	omsg << e.what() << std::endl \
  << "caught in " << __FILE__ << ":" << __LINE__ << std::endl ; \
  omsg << msg << std::endl; \
  throw Sundance::RuntimeError(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg)); \
}

#define SUNDANCE_BOUNDSCHECK(i, low, high, msg) \
{ \
  TEST_FOR_EXCEPTION( i < low || i > high, Sundance::RuntimeError, \
                     "Bounds violation: " << #i << "is out of range [" \
                      << #low << ", " << #high << "]") \
}

#define SUNDANCE_CHECK_ARRAY_SIZE_MATCH(a1, a2) \
  {\
    TEST_FOR_EXCEPTION(a1.size() != a2.size(), Sundance::RuntimeError, \
      "Mismatched array sizes: size(" << #a1 << ")=" << a1.size() \
      << " size(" << #a2 << ")=" << a2.size() << ". Expected equal sizes");\
  }



namespace Sundance
{
  /**
   * InternalError is thrown when an "impossible" case is detected
   * in Sundance's internals. Occurance of an InternalError indicates
   * either a bug in Sundance or runtime memory corruption that is
   * scrambling an object's virtual tables.
   */
  class InternalError : public std::logic_error
    {
    public:
      /** */
      InternalError(const std::string& msg);
    };

  /**
   * RuntimeError is an exception that occurs as a result of invalid
   * user-level code.
   */
  class RuntimeError : public std::runtime_error
    {
    public:
      /** */
      RuntimeError(const std::string& msg);
    };

  /**
   * BadSymbolicsError is thrown when a mathematically nonsensical
   * expression is detected. Unfortunately, it is possible to form expressions
   * that are legal C++ but illegal mathematics. For example, one can
   * happily divide a differential operator by another expression,
   * \code
   * Expr f = new UnknownFunction();
   * Expr dx = new Derivative(0);
   * Expr e = dx/f;
   * \endcode
   */
  class BadSymbolicsError : public RuntimeError
    {
    public:
      /** */
      BadSymbolicsError(const std::string& msg);
    };
}

                  
#endif  /* DOXYGEN_DEVELOPER_ONLY */  


#endif
