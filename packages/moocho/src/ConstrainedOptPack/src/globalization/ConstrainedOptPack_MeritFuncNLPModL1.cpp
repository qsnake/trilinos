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

#include "ConstrainedOptPack_MeritFuncNLPModL1.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "Teuchos_TestForException.hpp"

namespace ConstrainedOptPack {

MeritFuncNLPModL1::MeritFuncNLPModL1()
  : deriv_(0.0)
{}

// Overridden from MeritFuncNLP

value_type MeritFuncNLPModL1::value(
  value_type             f
  ,const Vector    *c
  ,const Vector    *h
  ,const Vector    *hl
  ,const Vector    *hu
  ) const
{
  TEST_FOR_EXCEPTION(
    h || hl || hu, std::logic_error
    ,"MeritFuncNLPModL1::value(...) : Error! general inequalities are not supported!" );
/*
  using DenseLinAlgPack::norm_1;
  return f + local_constr_term( mu_, c, "calc_deriv" );
*/
  TEST_FOR_EXCEPT(true); // ToDo: Write a reduction operator for the above operation
  return 0.0;
}

value_type MeritFuncNLPModL1::deriv() const
{
  return deriv_;
}

void MeritFuncNLPModL1::print_merit_func(
  std::ostream& out, const std::string& L
  ) const
{
  out
    << L << "*** Define a modified L1 merit funciton that uses different\n"
    << L << "*** penalty parameters for each constriant.\n"
    << L << "*** (assumes Gc_k'*d_k + c_k = 0):\n"
    << L << "phi(f,c) = f + sum( mu(j) * abs(c(j)), j = 1,...,m )\n"
    << L << "Dphi(x_k,d_k) = Gf_k' * d_k - sum( mu(j) * abs(c(j)), j = 1,...,m )\n";
}

// Overridden from MeritFuncNLPDirecDeriv

value_type MeritFuncNLPModL1::calc_deriv(
  const Vector    &Gf_k
  ,const Vector   *c_k
  ,const Vector   *h_k
  ,const Vector   *hl
  ,const Vector   *hu
  ,const Vector   &d_k
  )
{
  TEST_FOR_EXCEPTION(
    h_k || hl || hu, std::logic_error
    ,"MeritFuncNLPModL1::value(...) : Error! general inequalities are not supported!" );
/*
  using DenseLinAlgPack::dot; using DenseLinAlgPack::norm_1;
  return deriv_ = dot( Gf_k, d_k ) - local_constr_term( mu_, c_k, "calc_deriv" );
*/
  TEST_FOR_EXCEPT(true); // ToDo: Write a reduction operator for the above operation
  return 0.0;
}

// Overridden from MeritFuncPenaltyParam

void MeritFuncNLPModL1::set_space_c( const VectorSpace::space_ptr_t& space_c )
{
  mu_  = space_c->create_member();
  *mu_ = 0.0;
}

VectorMutable& MeritFuncNLPModL1::set_mu()
{
  return *mu_;
}

const Vector& MeritFuncNLPModL1::get_mu() const
{
  return *mu_;
}

}	// end namespace ConstrainedOptPack

/* ToDo: Write a reduction operator for the following!

namespace {

value_type local_constr_term( const DVector& mu, const DVectorSlice& c
  , const char func_name[] )
{
  if( mu.size() != c.size() ) {
    std::ostringstream omsg;
    omsg
      << "MeritFuncNLPModL1::" << func_name << "(...) : "
      << "Error, the sizes mu.size() == " << mu.size()
      << " != c.size() == " << c.size();
    throw ConstrainedOptPack::MeritFuncNLP::InvalidInitialization(omsg.str());
  }
  value_type r = 0.0;
  DVector::const_iterator
    mu_itr = mu.begin();
  DVectorSlice::const_iterator
    c_itr = c.begin();
  while( mu_itr != mu.end() )
    r += *mu_itr++ * ::fabs( *c_itr++ );
  return r;
}

}	// end namespace

*/
