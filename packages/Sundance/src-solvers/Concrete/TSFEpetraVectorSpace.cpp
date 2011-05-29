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

#include "TSFEpetraVectorSpace.hpp"
#include "TSFEpetraVector.hpp"
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

EpetraVectorSpace::EpetraVectorSpace(const RCP<const Epetra_Map>& m)
  : ScalarProdVectorSpaceBase<double>(),
    SpmdVectorSpaceBase<double>(),
    smallVecSpcFactory_(rcp(new DefaultSpmdVectorSpaceFactory<double>())),
    epetraMap_(m),
    comm_(epetraCommToTeuchosComm(m->Comm())),
    localSubDim_(epetraMap_->NumMyElements()),
    localOffset_(epetraMap_->MinMyGID())
{}


OrdType EpetraVectorSpace::dim() const 
{
  return epetraMap_->NumGlobalElements();
}

bool EpetraVectorSpace::isCompatible(const VectorSpaceBase<double>& other) const 
{
  const EpetraVectorSpace* epvs = dynamic_cast<const EpetraVectorSpace*>(&other);
  if (epvs != 0)
  {
    return epetraMap_->SameAs(*(epvs->epetraMap_));
  }
  return false;
}

RCP<const VectorSpaceFactoryBase<double> > 
EpetraVectorSpace::smallVecSpcFcty() const 
{
  return smallVecSpcFactory_;
}



// Overridden from VectorSpace

Teuchos::RCP<VectorBase<double> >
EpetraVectorSpace::createMember() const
{
//  cout << "creating vector" << std::endl;
  return rcp(new EpetraVector(rcp(this, false)));
}



Teuchos::RCP<MultiVectorBase<double> >
EpetraVectorSpace::createMembers(int n) const
{
  RCP<const VectorSpaceBase<double> > self = rcp(this, false);
  RCP<const VectorSpaceBase<double> > small 
    = Thyra::defaultSpmdVectorSpace<double>(n);
  Array<RCP<VectorBase<double> > > vecs(n);
  for (int i=0; i<vecs.size(); i++)
    {
      vecs[i] = createMember();
    }
  return rcp(
    new Thyra::DefaultColumnwiseMultiVector<double>(
      self, small,
#ifndef TRILINOS_8
      vecs
#else
      &(vecs[0])
#endif
      )
    );
}



Teuchos::RCP< const VectorSpaceBase<double> >
EpetraVectorSpace::clone() const
{
  return Teuchos::rcp(new EpetraVectorSpace(epetraMap_));
}



string EpetraVectorSpace::description() const
{
  std::string rtn = "EpetraVS[d=" + Teuchos::toString(dim());
  if (localSubDim() != dim()) rtn += ", local="
    + Teuchos::toString(localSubDim());
  rtn += "]";
  return rtn;
}



Teuchos::RCP<const Teuchos::Comm<OrdType> > 
EpetraVectorSpace::epetraCommToTeuchosComm(const Epetra_Comm& epComm) const 
{
  RCP<const Comm<OrdType> > rtn;

#ifdef HAVE_MPI
  const Epetra_MpiComm* mpiComm 
    = dynamic_cast<const Epetra_MpiComm*>(&epComm);
#endif

  const Epetra_SerialComm* serialComm 
    = dynamic_cast<const Epetra_SerialComm*>(&epComm);

  if (serialComm != 0)
  {
    rtn  = rcp(new SerialComm<OrdType>());
  }
#ifdef HAVE_MPI
  else if (mpiComm != 0)
  {
    MPI_Comm rawMpiComm = mpiComm->GetMpiComm();
    RCP<const OpaqueWrapper<MPI_Comm> > ptr 
      = rcp(new OpaqueWrapper<MPI_Comm>(rawMpiComm));
    rtn  = rcp(new MpiComm<OrdType>(ptr));
  }
#endif
  else
  {
    TEST_FOR_EXCEPTION(true, std::runtime_error, "Epetra_Comm is neither "
      "a SerialComm or MpiComm");
  }
  return rtn;
}




