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

#ifndef SET_OPTIONS_FROM_STREAM_NODE_H
#define SET_OPTIONS_FROM_STREAM_NODE_H

#include "OptionsFromStreamPack_SetOptionsFromStream.hpp"
#include "OptionsFromStreamPack_StringToIntMap.hpp"

namespace OptionsFromStreamPack {

/** \brief Node class for setting options from a stream.
  *
  * This class uses the template method pattern to
  * delegate the setting of options.
  */
class SetOptionsFromStreamNode: public SetOptionsFromStream {
public:

  /** \brief Constructs with the name of the options group and the names
    * of the options.
    *
    *	@param	options_group	The name of the options group to access
    *	@param	num_options		The number of options in the opitons
    *							group.
    *	@param	option_name		An array (length num_options) containing
    *							the names of the options.
    *	@param	exists_optional	Specifies if the options group must exist.
    */
  SetOptionsFromStreamNode( const std::string& options_group
    , int num_options, const char* option_names[]
    , bool exists_optional = true );

  /** \brief Overridden from SetOptionsFromStream and calls setOption(...).
    *
    * The options group #options_group# is used.  If this options
    * group does not exist and #exists_optional# == false then
    * an #std::invalid_argument# exception will be thrown.
    */
  void set_options( const OptionsFromStream& options );

protected:

  /** \brief To be overridden by the subclass to set an option given
    * its integer position and the option value.
    *
    * The integer possition returned is the possition of the option
    * in option_names[option_num] that was passed to the constructor.
    */
  virtual void setOption( int option_num, const std::string& option_value ) = 0;

private:
  StringToIntMap	name_map_;
  bool			exists_optional_;

};	// end class SetOptionsFromStreamNode

}	// end namespace OptionsFromStreamPack

#endif	// SET_OPTIONS_FROM_STREAM_NODE_H
