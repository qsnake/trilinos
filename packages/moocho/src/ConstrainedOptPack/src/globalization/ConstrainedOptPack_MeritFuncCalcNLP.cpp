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

#include "ConstrainedOptPack_MeritFuncCalcNLP.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"

namespace ConstrainedOptPack {

MeritFuncCalcNLP::MeritFuncCalcNLP( const MeritFuncNLP* phi, const NLP* nlp )
  : phi_(phi), nlp_(nlp)
{}

value_type MeritFuncCalcNLP::operator()(const Vector& x) const
{
  const size_type
    m  = nlp().m(),
    ns = nlp().ns();
  nlp().calc_f(x);
  if(m)  nlp().calc_c(x,false);
  return phi().value(
    nlp().f()
    ,m  ? &nlp().c()  : NULL
    ,NULL  // h
    ,NULL  // hl
    ,NULL  // hu
    );
/* RAB: 20020112: ToDo: Get this working
  if(m)  nlp().calc_c_breve(x,false);
  if(ns) nlp().calc_h_breve(x,false);
  return phi().value(
    nlp().f()
    ,m  ? &nlp().c_breve()  : NULL
    ,ns ? &nlp().h_breve()  : NULL
    ,ns ? &nlp().hl_breve() : NULL
    ,ns ? &nlp().hu_breve() : NULL
    );
*/
}

value_type MeritFuncCalcNLP::deriv() const {
  return phi().deriv();
}

void MeritFuncCalcNLP::print_merit_func(
  std::ostream& out, const std::string& L
  ) const
{
  out	<< L << "*** MeritFuncCalcNLP\n"
    << L << "f = f(x), c = c_breve(x_breve), h = h_breve(x_breve)\n";
  phi().print_merit_func(out,L);
}

}	// end namespace ConstrainedOptPack
