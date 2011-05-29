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

#include "TSFSerialVectorSpace.hpp"
#include "TSFSerialVector.hpp"
#include "Teuchos_Utils.hpp"

#include "Teuchos_DefaultSerialComm.hpp"
#include "Thyra_DefaultSpmdVectorSpaceFactory.hpp"
#include "Thyra_DefaultSpmdVectorSpace_decl.hpp"
#include "Thyra_DefaultColumnwiseMultiVector.hpp"
#include "SundanceOut.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_Comm.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "Teuchos_DefaultMpiComm.hpp"
#endif


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorSpaceImpl.hpp"
#endif



using namespace TSFExtended;
using namespace Thyra;

using Teuchos::RCP;

SerialVectorSpace::SerialVectorSpace(int dim)
  : DefaultSpmdVectorSpace<double>(dim),
    smallVecSpcFactory_(rcp(new DefaultSpmdVectorSpaceFactory<double>()))
{
}


RCP<const VectorSpaceFactoryBase<double> > 
SerialVectorSpace::smallVecSpcFcty() const 
{
  return smallVecSpcFactory_;
}

bool SerialVectorSpace::isCompatible(const VectorSpaceBase<double>& other) const 
{
  const SerialVectorSpace* rvs 
    = dynamic_cast<const SerialVectorSpace*>(&other);
  if (rvs != 0)
  {
    return this->dim() == rvs->dim();
  }
  return false;
}


// Overridden from VectorSpace

Teuchos::RCP<VectorBase<double> >
SerialVectorSpace::createMember() const
{
  return rcp(new SerialVector(rcp(this, false)));
}



Teuchos::RCP<MultiVectorBase<double> >
SerialVectorSpace::createMembers(int n) const
{
  RCP<const VectorSpaceBase<double> > self = rcp(this, false);
  RCP<const VectorSpaceBase<double> > small 
    = rcp(new SerialVectorSpace(n));
  Array<RCP<VectorBase<double> > > vecs(n);
  for (int i=0; i<vecs.size(); i++)
  {
    vecs[i] = createMember();
  }
  return rcp(
    new Thyra::DefaultColumnwiseMultiVector<double>(
      self, small,
      vecs
      )
    );
}



Teuchos::RCP< const VectorSpaceBase<double> >
SerialVectorSpace::clone() const
{
  return Teuchos::rcp(new SerialVectorSpace(this->dim()));
}



string SerialVectorSpace::description() const
{
  std::string rtn = "SerialVS[d=" + Teuchos::toString(this->dim()) + "]";
  return rtn;
}




