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

#include "NLPInterfacePack_NLPObjGrad.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"

namespace {
  const char name_Gf[] = "Gf";
} // end namespace

namespace NLPInterfacePack {

// constructors

NLPObjGrad::NLPObjGrad()
  : Gf_(NULL)
{}

void NLPObjGrad::initialize(bool test_setup) {
  num_Gf_evals_ = 0;
  NLP::initialize(test_setup);
}

// Information

bool NLPObjGrad::supports_Gf() const
{
  return true;
}

bool NLPObjGrad::supports_Gf_prod() const
{
  return false;
}

// <<std aggr>> members for Gf

void NLPObjGrad::set_Gf(VectorMutable* Gf)
{
  Gf_ = Gf;
}

AbstractLinAlgPack::VectorMutable* NLPObjGrad::get_Gf()
{
  return StandardCompositionRelationshipsPack::get_role_name(Gf_, false, name_Gf);
}

AbstractLinAlgPack::VectorMutable& NLPObjGrad::Gf()
{
  return StandardCompositionRelationshipsPack::role_name(Gf_, false, name_Gf);
}

const AbstractLinAlgPack::Vector& NLPObjGrad::Gf() const
{
  return StandardCompositionRelationshipsPack::role_name(Gf_, false, name_Gf);
}

void NLPObjGrad::unset_quantities()
{
  NLP::unset_quantities();
  Gf_ = NULL;
}

// calculations

void NLPObjGrad::calc_Gf(const Vector& x, bool newx) const
{
  StandardCompositionRelationshipsPack::assert_role_name_set(Gf_, "NLP::calc_Gf()", name_Gf);
  imp_calc_Gf(x,newx,obj_grad_info());
  num_Gf_evals_++;
}

value_type NLPObjGrad::calc_Gf_prod(const Vector& x, const Vector& d, bool newx) const
{
  TEST_FOR_EXCEPTION(
    true,std::logic_error
    ,"Error, the function calc_Gf_prod(...) is not implemented for the class "
    << typeName(*this) << "!"
    );

  //execution should never reach this point, but compilers expect a non-void
  //function to return something. So we'll create a dummy value to use in a
  //return statement.
  //(a better design would not require function bodies for unimplemented
  //functions like this...)
  value_type* dummy = NULL;
  return(*dummy);
}

size_type NLPObjGrad::num_Gf_evals() const
{
  return num_Gf_evals_;
}

}	// end namespace NLPInterfacePack 
