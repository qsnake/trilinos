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

#include <math.h>
#include <iostream>
#include <limits>

#include "NLPInterfacePack_NLPBarrier.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "Teuchos_TestForException.hpp"

namespace NLPInterfacePack {

NLPBarrier::NLPBarrier()
  :
  barrier_term_(0.0),
  objective_term_(0.0),
  nlp_(Teuchos::null)
  {
  }


void NLPBarrier::InitializeFromNLP(
  Teuchos::RCP<NLP> original_nlp
  )
  {
  TEST_FOR_EXCEPTION(
    !original_nlp.get(),
    std::logic_error,
    "null nlp passed to NLPBarrier decorator"
    );

  nlp_ = Teuchos::rcp_dynamic_cast<NLPObjGrad>(original_nlp);

  TEST_FOR_EXCEPTION(
    !nlp_.get(),
    std::logic_error,
    "non NLPObjGrad NLP passed to NLPBarrier decorator"
    );
  }

void NLPBarrier::mu(const value_type mu)
  {
  mu_ = mu;
  }

value_type NLPBarrier::barrier_term() const
  {
  return barrier_term_;
  }

value_type NLPBarrier::objective_term() const
  {
  return objective_term_;
  }

const Teuchos::RCP<Vector> NLPBarrier::grad_barrier_term() const
  {
  return grad_barrier_term_;
  }

const Teuchos::RCP<Vector>  NLPBarrier::grad_objective_term() const
  {
  return grad_objective_term_;
  }


void NLPBarrier::calc_f(const Vector& x, bool newx) const
  {
  nlp_->calc_f(x, newx);
  value_type* f_val = nlp_->get_f();

  objective_term_ = *f_val;
  barrier_term_   = CalculateBarrierTerm(x);

  (*f_val) += barrier_term_;
  }

void NLPBarrier::calc_Gf(const Vector& x, bool newx) const
  {
  using AbstractLinAlgPack::inv_of_difference;

     nlp_->calc_Gf(x, newx);
  grad_objective_term_ = nlp_->get_Gf()->clone();

  //std::cout << "grad_objective_term=\n";
  //grad_objective_term_->output(std::cout);

  if (!grad_barrier_term_temp_.get())
    { grad_barrier_term_temp_ = grad_objective_term_->space().create_member(); }

  if (!grad_barrier_term_.get())
    { grad_barrier_term_ = grad_objective_term_->space().create_member(); }	

  *grad_barrier_term_temp_ = 0.0;
  *grad_barrier_term_ = 0.0;	

  inv_of_difference(mu_, nlp_->xu(), x, grad_barrier_term_.get());
   //std::cout << "mu*invXU=\n";
  //grad_barrier_term_->output(std::cout);

  inv_of_difference(mu_, x, nlp_->xl(), grad_barrier_term_temp_.get());
   //std::cout << "mu*invXL=\n";
  //grad_barrier_term_temp_->output(std::cout);

  grad_barrier_term_->axpy(-1.0, *grad_barrier_term_temp_);

  nlp_->get_Gf()->axpy(1.0, *grad_barrier_term_);

  //std::cout << "grad_objective_term with barrier=\n";
  //nlp_->get_Gf()->output(std::cout);
  }

void NLPBarrier::imp_calc_f(
  const Vector& x, 
  bool newx, 
  const ZeroOrderInfo& zero_order_info
  ) const
  {
  TEST_FOR_EXCEPT( !( false && !"This should never get called." ) );
  }

void NLPBarrier::imp_calc_c(
  const Vector& x, 
  bool newx, 
  const ZeroOrderInfo& zero_order_info
  ) const
  {
  TEST_FOR_EXCEPT( !( false && !"This should never get called." ) );
  }

void NLPBarrier::imp_calc_c_breve(
  const Vector& x, 
  bool newx, 
  const ZeroOrderInfo& zero_order_info_breve
  ) const
  {	
  TEST_FOR_EXCEPT( !( false && !"This should never get called." ) );
  }

void NLPBarrier::imp_calc_h_breve(
  const Vector& x, 
  bool newx, 
  const ZeroOrderInfo& zero_order_info_breve
  ) const
  {	
  TEST_FOR_EXCEPT( !( false && !"This should never get called." ) );
  }

void NLPBarrier::imp_calc_Gf(
  const Vector& x,
  bool newx, 
  const ObjGradInfo& obj_grad_info
  ) const
  {
  TEST_FOR_EXCEPT( !( false && !"This should never get called." ) );
  }


value_type NLPBarrier::CalculateBarrierTerm(const Vector& x) const
  {
  using AbstractLinAlgPack::log_bound_barrier;
  barrier_term_ = log_bound_barrier(x, xl(), xu());
//	std::cerr << "NLPBarrier::CalculateBarrierTerm(x) : (1) barrier_term_ = " << barrier_term_ << std::endl;
  barrier_term_ *= -mu_;
//	std::cerr << "NLPBarrier::CalculateBarrierTerm(x) : (2) barrier_term_ = " << barrier_term_ << std::endl;
  return barrier_term_;
  }

}	// end namespace NLPInterfacePack
