// @HEADER
// ***********************************************************************
// 
//         Optika: A Tool For Developing Parameter Obtaining GUIs
//                Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, with Sandia Corporation, the 
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Kurtis Nusbaum (klnusbaum@gmail.com) 
// 
// ***********************************************************************
// @HEADER
#ifndef OPTIKA_INVALIDCONDITIONEXCEPTION_HPP_
#define OPTIKA_INVALIDCONDITIONEXCEPTION_HPP_
#include <stdexcept>

namespace Optika {

/**
 * Thrown when some aspect of a Condition has been determined to be invalid.
 */
class InvalidConditionException : public std::logic_error{
public: 
	/**
	 * Constructs an InvalidConditionException
	 *
	 * @param what_arg The error message to be associated with this error.
	 */
	InvalidConditionException(const std::string& what_arg):std::logic_error(what_arg){}
};

}
#endif //OPTIKA_INVALIDCONDITIONEXCEPTION_HPP_
