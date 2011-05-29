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

#include "TSFDenseSerialMatrixFactory.hpp"
#include "TSFDenseSerialMatrix.hpp"
#include "TSFSerialVector.hpp"
#include "TSFVectorSpaceDecl.hpp"  
#include "TSFVectorDecl.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_MPIComm.hpp"
#include "TSFLinearOperatorDecl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#include "TSFVectorImpl.hpp"
#endif


using namespace TSFExtended;
using namespace Teuchos;
using namespace Thyra;

DenseSerialMatrixFactory::DenseSerialMatrixFactory(
  const RCP<const SerialVectorSpace>& domain,
  const RCP<const SerialVectorSpace>& range)
  : range_(range),
    domain_(domain)
{}


LinearOperator<double> DenseSerialMatrixFactory::createMatrix() const
{
  RCP<LinearOpBase<double> > A 
    = rcp(new DenseSerialMatrix(domain_, range_));
  return A;
}
