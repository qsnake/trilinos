/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// **********************************************************************/
 /* @HEADER@ */

#ifndef SUNDANCE_NAMED_OBJECT_HPP
#define SUNDANCE_NAMED_OBJECT_HPP

#include "SundanceDefs.hpp"
#include <string>

namespace Sundance
{
/**
 * NamedObject provides an interface for getting and setting 
 * names of an object.
 *
 * @author Kevin Long (kevin.long@ttu.edu)
 */
class NamedObject 
{
public:
  /** */
  NamedObject() : name_("Anonymous[]") {}

  /** */
  virtual ~NamedObject() {}
      
  /** */
  void setName(const std::string& name) {name_ = name;}

  /** */
  std::string name() const 
    {return name_;}

private:
  std::string name_;
      
};


}


#endif
