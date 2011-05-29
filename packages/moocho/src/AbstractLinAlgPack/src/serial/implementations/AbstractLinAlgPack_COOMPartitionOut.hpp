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
// Output stream operators for Partition<> and TransposedPartition<>

#ifndef COOM_PARTITION_OUT_H
#define COOM_PARTITION_OUT_H

#include "AbstractLinAlgPack_COOMatrixTmplOutFunc.hpp"

namespace AbstractLinAlgPack {

/** \brief Partition<> output stream operator.
  *
  * This operator function calls the function output_COOM(os,part,0).
  */
template <class T_Indice, class T_Value>
inline std::ostream& operator<<(std::ostream& os
  , const COOMatrixPartitionedViewUtilityPack::Partition<T_Indice,T_Value>& part) {
  return output_COOM(os,part,0);
}

/** \brief TransposedPartition<> output stream operator.
  *
  * This operator function calls the function output_COOM(os,trans_part,0).
  */
template <class T_Indice, class T_Value>
inline std::ostream& operator<<(std::ostream& os
  , const COOMatrixPartitionedViewUtilityPack::TransposedPartition<T_Indice,T_Value>& trans_part)
{
  return output_COOM(os,trans_part,0);	
}

}	// end namespace AbstractLinAlgPack

#endif // VECTOROUT_H
