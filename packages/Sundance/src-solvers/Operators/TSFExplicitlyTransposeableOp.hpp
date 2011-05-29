/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// **********************************************************************/
/* @HEADER@ */

#ifndef TSFEXPLICITLYTRANSPOSEABLEOP_HPP
#define TSFEXPLICITLYTRANSPOSEABLEOP_HPP

#include "SundanceDefs.hpp"
#include "TSFLinearOperatorDecl.hpp"


namespace TSFExtended
{

  /** 
   * Base interface for operators whose transpose can be formed
   * explicitly. Note: if a transpose is formed, it will be an independent
   * object in the sense that any changes made to either the original object
   * or the transpose are not propagated to the other.
   */
  template <class Scalar>
  class ExplicitlyTransposeableOp 
    {
    public:
      /** Virtual dtor */
      virtual ~ExplicitlyTransposeableOp(){;}

      /** Form a new object that is the transpose of this operator */
      virtual LinearOperator<Scalar> formTranspose() const = 0 ;
      
    private:
    };
}


#endif
