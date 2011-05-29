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

#ifndef CALC_FINITE_DIFF_PROD_SET_OPTIONS_H
#define CALC_FINITE_DIFF_PROD_SET_OPTIONS_H

#include "NLPInterfacePack_CalcFiniteDiffProd.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace NLPInterfacePack {

/** \brief Set options for \c CalcFiniteDiffProd from an
 * \c OptionsFromStream object.
 *
 * The default options group name is CalcFiniteDiffProd.
 *
 * The options group is:
 *
 \verbatim
 options_group CalcFiniteDiffProdSetOptions {
 *    fd_method_order = FD_ORDER_ONE;
 *    fd_method_order = FD_ORDER_TWO;
 *    fd_method_order = FD_ORDER_TWO_CENTRAL;
 *    fd_method_order = FD_ORDER_TWO_AUTO;
 *    fd_method_order = FD_ORDER_FOUR;
 *    fd_method_order = FD_ORDER_FOUR_CENTRAL;
 *    fd_method_order = FD_ORDER_FOUR_AUTO; *** (Default)
 *    fd_step_select = FD_STEP_ABSOLUTE; *** (Default)
 *    fd_step_select = FD_STEP_RELATIVE;
 *    fd_step_size = -1.0; *** (default)
 *    fd_step_size_min = -1.0; *** (default)
 *    fd_step_size_f = -1.0; *** (default)
 *    fd_step_size_c = -1.0; *** (default)
 }
 \endverbatim
 */
class CalcFiniteDiffProdSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      CalcFiniteDiffProd >
{
public:

  /** \brief . */
  CalcFiniteDiffProdSetOptions(
    CalcFiniteDiffProd* target = 0
    ,const char opt_grp_name[] = "CalcFiniteDiffProd" );
  
protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class CalcFiniteDiffProdSetOptions

}	// end namespace NLPInterfacePack

#endif	// CALC_FINITE_DIFF_PROD_SET_OPTIONS_H
