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

#ifndef QUASI_NEWTON_STATS_HH
#define QUASI_NEWTON_STATS_HH

#include "MoochoPack_QuasiNewtonStats.hpp"
#include "IterationPack_CastIQMember.hpp"

namespace MoochoPack {

/// Name given to the quasi-Newton updating staistics iteration quantity
extern const std::string quasi_newton_stats_name;

/** \brief Class for object that attempts to return an IterQuantityAccess<QuasiNewtonStats>
  * from an AlgorithmState object with the name quasi_newton_stats_name.
  */
class quasi_newton_stats_iq_member
  : public CastIQMember<QuasiNewtonStats>
{
public:
    quasi_newton_stats_iq_member()
      : CastIQMember<QuasiNewtonStats>(quasi_newton_stats_name)
    {}
};

}	// end namespace MoochoPack

#endif	// QUASI_NEWTON_STATS_HH
