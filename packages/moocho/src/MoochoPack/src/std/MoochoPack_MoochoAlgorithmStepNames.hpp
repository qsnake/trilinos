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

#ifndef RSQP_ALGORITHM_STEP_NAMES_H
#define RSQP_ALGORITHM_STEP_NAMES_H

#include <string>

#include "MoochoPack_Types.hpp"

namespace MoochoPack {

/** @name Names for MOOCHOsteps */
//@{

const std::string EvalNewPoint_name                 = "EvalNewPoint";
const std::string ReducedGradient_name              = "ReducedGradient";
const std::string ReducedHessian_name               = "ReducedHessian";
const std::string QuasiNormalStep_name              = "QuasiNormalStep";
const std::string TangentialStep_name               = "TangentialStep";
const std::string SearchDirec_name                  = "SearchDirec";
const std::string LineSearch_name                   = "LineSearch";
const std::string CheckConvergence_name             = "CheckConvergence";

const std::string CalcLambdaIndep_name              = "CalcLambdaIndep";
const std::string CalcReducedGradLagrangian_name    = "CalcReducedGradLagrangian";
const std::string CheckSkipBFGSUpdate_name          = "CheckSkipBFGSUpdate";
const std::string CalcDFromYPYZPZ_name				= "CalcDFromYPYZPZ";

//@}
}	// end namespace MoochoPack

#endif // RSQP_ALGORITHM_STEP_NAMES_H
