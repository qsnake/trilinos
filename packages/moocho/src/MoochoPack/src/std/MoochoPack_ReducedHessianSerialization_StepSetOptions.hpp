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

#ifndef REDUCED_HESSIAN_SERIALIZATION_STEP_SET_OPTIONS_H
#define REDUCED_HESSIAN_SERIALIZATION_STEP_SET_OPTIONS_H

#include "MoochoPack_ReducedHessianSerialization_Step.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Set options for ReducedHessianSerialization_Step from an \c OptionsFromStream object.
 *
 * The default options group name is EvalNewPointStd.
 *
 * The options group is:
 *
 \verbatim
    options_group ReducedHessianSerialization {
      reduced_hessian_input_file_name   = "reduced_hessian.in";
      reduced_hessian_output_file_name  = "reduced_hessian.out";
    }
 \verbatim
 *
 * See the file Moocho.opt.NLPAlgoConfigMamaJama
 */
class ReducedHessianSerialization_StepSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      ReducedHessianSerialization_Step >
{
public:

  /** \brief . */
  ReducedHessianSerialization_StepSetOptions(
      ReducedHessianSerialization_Step* target = 0
    , const char opt_grp_name[] = "ReducedHessianSerialization" );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class ReducedHessianSerialization_StepSetOptions

}	// end namespace MoochoPack

#endif	// REDUCED_HESSIAN_SERIALIZATION_STEP_SET_OPTIONS_H
