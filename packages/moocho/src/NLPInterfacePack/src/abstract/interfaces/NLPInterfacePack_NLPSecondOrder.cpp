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

#include "NLPInterfacePack_NLPSecondOrder.hpp"
#include "Teuchos_TestForException.hpp"

namespace {
  const char name_HL[] = "HL";
}

namespace NLPInterfacePack {

// constructors

NLPSecondOrder::NLPSecondOrder()
  : HL_(NULL)
{}


void NLPSecondOrder::initialize(bool test_setup) {
  num_HL_evals_ = 0;
  NLPFirstOrder::initialize(test_setup);
}

// <<std aggr>> members for HL

void NLPSecondOrder::set_HL(MatrixSymOp* HL)
{
  HL_ = HL;
}

MatrixSymOp* NLPSecondOrder::get_HL()
{
  return StandardCompositionRelationshipsPack::get_role_name(HL_, false, name_HL);
}

MatrixSymOp& NLPSecondOrder::HL()
{
  return StandardCompositionRelationshipsPack::role_name(HL_, false, name_HL);
}

const MatrixSymOp& NLPSecondOrder::HL() const
{
  return StandardCompositionRelationshipsPack::role_name(HL_, false, name_HL);
}

void NLPSecondOrder::unset_quantities()
{
  NLPFirstOrder::unset_quantities();
  HL_ = NULL;
}

// calculations

void NLPSecondOrder::calc_HL(
  const Vector& x, const Vector* lambda, bool newpoint
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION( lambda  && this->m()  == 0, std::logic_error, "" );
#endif
  StandardCompositionRelationshipsPack::assert_role_name_set(HL_, "NLP::calc_HL()", name_HL);
  imp_calc_HL(x,lambda,newpoint,second_order_info());
}

size_type NLPSecondOrder::num_HL_evals() const
{
  return num_HL_evals_;
}

} // namespace NLPInterfacePack
