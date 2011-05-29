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

#include <ostream>

#include "IterationPack_TestIterationPack.hpp"
#include "TestingHelperPack_update_success.hpp"

bool IterationPack::TestingPack::TestIterationPack(std::ostream* out)
{
  using TestingHelperPack::update_success;
  using namespace IterationPack::TestingPack;

  bool success = true;
  update_success( TestIterQuantityAccessContiguous(out), &success );
  update_success( TestAlgorithmState(out), &success );
  update_success( TestAlgorithm(out), &success );

  if(out) {
    if(success) {
      *out << "\n*** Congradulations, IterationPack seems to"
          " check out.\n";
    }
    else {
      *out << "\n*** Oops, at least one of the above tests for IterationPack "
          " did not return the expected results\n";
    }
  } 

  return success;
}
