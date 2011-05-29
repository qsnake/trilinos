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

#ifndef ACT_SET_STATS_HH
#define ACT_SET_STATS_HH

#include "MoochoPack_ActSetStats.hpp"
#include "IterationPack_CastIQMember.hpp"

namespace MoochoPack {

/// Name given to the active set statistics iteration quantity
extern const std::string act_set_stats_name;

/** \brief Class for object that attempts to return an IterQuantityAccess<ActSetStats>
  * from an AlgorithmState object with the name act_set_stats_name.
  */
class act_set_stats_iq_member
  : public CastIQMember<ActSetStats>
{
public:
    act_set_stats_iq_member()
      : CastIQMember<ActSetStats>(act_set_stats_name)
    {}
};

}	// end namespace MoochoPack

#endif	// ACT_SET_STATS_HH
