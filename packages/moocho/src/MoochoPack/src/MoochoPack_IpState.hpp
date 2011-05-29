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
//

#if !defined IP_STATE_H
#define IP_STATE_H

#include "MoochoPack_NLPAlgoState.hpp"

namespace MoochoPack {

// Iteration Quantity Strings
extern const std::string barrier_parameter_name;
extern const std::string barrier_obj_name;
extern const std::string grad_barrier_obj_name;
extern const std::string e_tol_name;
extern const std::string comp_err_mu_name;
extern const std::string Vu_name;
extern const std::string Vl_name;
extern const std::string invXu_name;
extern const std::string invXl_name;
extern const std::string rHB_name;
extern const std::string B_name;
extern const std::string Sigma_name;
extern const std::string w_sigma_name;
extern const std::string dvl_name;
extern const std::string dvu_name;
extern const std::string alpha_vl_name;
extern const std::string alpha_vu_name;


class IpState 
  : public MoochoPack::NLPAlgoState
  {

  public:
    ///********** Iteration Quantities **************

    /// mu: barrier parameter
    STATE_SCALAR_IQ_DECL(barrier_parameter)

    /// barrier_obj: objective value with 
    //   barrier term included
    STATE_SCALAR_IQ_DECL(barrier_obj)

    /// grad_barrier_obj: gradient of the objective
    //   with barrier term included
    STATE_VECTOR_IQ_DECL(grad_barrier_obj)

    /// e_tol: current error tolerance for inner loop
    STATE_SCALAR_IQ_DECL(e_tol)

    /// comp_err_mu: perturbed complementarity error for barrier sub problem
    STATE_SCALAR_IQ_DECL(comp_err_mu)

    /// Vu - diagonal matrix of upper bound multipliers
    STATE_IQ_DECL(MatrixSymDiagStd, Vu)

    /// Vl - diagonal matrix of lower bound multipliers
    STATE_IQ_DECL(MatrixSymDiagStd, Vl)

    /// invXu - (Xu)^-1 - matrix of 1/(xu-x) diagonal
    STATE_IQ_DECL(MatrixSymDiagStd, invXu)

    /// invXl - (Xl)^-1 - matrix of 1/(x-xl) diagonal
    STATE_IQ_DECL(MatrixSymDiagStd, invXl)

    /// rHB - reduced Hessian of the barrier term (Z_Sigma_Z)
    STATE_IQ_DECL(MatrixSymOp, rHB)

    /// B - overall reduced 'Hessian' (Z_W_Z+Z_Sigma_Z)
    STATE_IQ_DECL(MatrixSymOp, B)

    /// Full space Sigma (invXl*Vl-invXu*Vu)
    STATE_IQ_DECL(MatrixSymDiagStd, Sigma)

    /// w_sigma:  crossterm correction for sigma (Z' * Sigma * Y * py)
    STATE_VECTOR_IQ_DECL(w_sigma) 

    /// dvl:  Search direction for lower bound multipliers ( n x 1 )
    STATE_VECTOR_IQ_DECL(dvl)

    /// dvu:  Search direction for upper bound multipliers ( n x 1 )
    STATE_VECTOR_IQ_DECL(dvu)

    /// alpha_vl: step size for vl
    STATE_SCALAR_IQ_DECL(alpha_vl)

    /// alpha_vl: step size for vu
    STATE_SCALAR_IQ_DECL(alpha_vu)

    /** \brief Construct
     *
     * 
     */
    IpState(
      const decomp_sys_ptr_t& decomp_sys   = Teuchos::null
      ,const vec_space_ptr_t& space_x      = Teuchos::null
      ,const vec_space_ptr_t& space_c      = Teuchos::null
      ,const vec_space_ptr_t& space_range  = Teuchos::null
      ,const vec_space_ptr_t& space_null   = Teuchos::null
      );

    virtual ~IpState();

  }; // end class IpState

} // end namespace MoochoPack



#endif // if !defined IP_STATE_H
