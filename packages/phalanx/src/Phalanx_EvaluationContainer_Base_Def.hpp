// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#ifndef PHX_SCALAR_CONTAINER_BASE_DEF_HPP
#define PHX_SCALAR_CONTAINER_BASE_DEF_HPP

#include "Teuchos_TestForException.hpp"
// **************************************************************************
template<typename Traits>
PHX::EvaluationContainerBase<Traits>::EvaluationContainerBase()
{

}

// **************************************************************************
template<typename Traits>
PHX::EvaluationContainerBase<Traits>::~EvaluationContainerBase()
{

}

// **************************************************************************
template<typename Traits>
void PHX::EvaluationContainerBase<Traits>::
requireField(const PHX::FieldTag& f) 
{ 
  vp_manager_.requireField(f);
}

// **************************************************************************
template<typename Traits>
void PHX::EvaluationContainerBase<Traits>::
registerEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits> >& e) 
{ 
  vp_manager_.registerEvaluator(e);
}

// **************************************************************************
template<typename Traits>
void PHX::EvaluationContainerBase<Traits>::
writeGraphvizFile(const std::string filename,
		  bool writeEvaluatedFields,
		  bool writeDependentFields,
		  bool debugRegisteredEvaluators) const 
{ 
  vp_manager_.writeGraphvizFile(filename, 
				writeEvaluatedFields, 
				writeDependentFields, 
				debugRegisteredEvaluators);
}
    
// **************************************************************************
template<typename Traits>
std::ostream&
PHX::operator<<(std::ostream& os, const PHX::EvaluationContainerBase<Traits>& sc)
{ 
  sc.print(os);
  return os;
}

// **************************************************************************

#endif
