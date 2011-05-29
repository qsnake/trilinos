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

#include "TSFVectorSpace2EpetraMap.hpp"
#include "TSFEpetraVectorSpace.hpp"
#include "TSFEpetraVector.hpp"
#include "Teuchos_Utils.hpp"
#include "Teuchos_DefaultSerialComm.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorSpaceImpl.hpp"
#endif

#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_SerialComm.h"


#ifdef TRILINOS_6
#include "Thyra_MultiVectorCols.hpp"
#ifdef HAVE_MPI
#include "Thyra_MPIVectorSpaceBase.hpp"
#endif
#include "Thyra_SerialVectorSpaceStd.hpp"
#define DefaultSerialVectorSpace SerialVectorSpaceStd
#define DefaultColumnwiseMultiVector MultiVectorCols
#else
#include "Thyra_DefaultSpmdVector.hpp"
#include "Thyra_SpmdVectorSpaceDefaultBase.hpp"
#include "Thyra_DefaultColumnwiseMultiVector.hpp"
#define MPIVectorSpaceBase SpmdVectorSpaceDefaultBase
#endif

namespace TSFExtended {
  using namespace Teuchos;
  using namespace Thyra;


  RCP<const Epetra_Map> tsfVectorSpace2EpetraMap(const VectorSpace<double>& tsfSpace)
  {
    const EpetraVectorSpace* ep 
      = dynamic_cast<const EpetraVectorSpace*>(tsfSpace.ptr().get());
    // See if we have an honest EpetraVectorSpace
    if (ep) 
      {
	return ep->epetraMap();
      }
    // Otherwise, make a EpetraMap with the same structure
    Array<int> globIndices;
    for (SequentialIterator<double> i = tsfSpace.begin(); 
	 i != tsfSpace.end(); i++)
      {
	globIndices.append(i.globalIndex());
      }
    int dim = tsfSpace.dim();

    RCP<Epetra_Comm> comm;
    TSFExtended::getComm(tsfSpace, comm);

    RCP<const Epetra_Map> rtn = rcp(new Epetra_Map(dim, globIndices.size(),
						     &(globIndices[0]),
						     0, *comm));
    return rtn;
  }


void getComm(const TSFExtended::VectorSpace<double>& tsfSpace,
    Teuchos::RCP<Epetra_Comm>& comm)
  {
#ifdef HAVE_MPI
    if (tsfSpace.numBlocks()==1)
      {
	const MPIVectorSpaceBase<double>* mv 
	  = dynamic_cast<const MPIVectorSpaceBase<double>*>(tsfSpace.getBlock(0).ptr().get());
	RCP<const Teuchos::Comm<OrdType> > tc = mv->getComm();
	const Teuchos::MpiComm<int>* mc 
	  = dynamic_cast<const Teuchos::MpiComm<int>*>(tc.get());
	const Teuchos::SerialComm<int>* sc 
	  = dynamic_cast<const Teuchos::SerialComm<int>*>(tc.get());
	if (mc != 0)
	  {
	    comm = rcp(new Epetra_MpiComm(*(mc->getRawMpiComm())));
	  }
	else if (sc != 0) 
	  {
	    comm = rcp(new Epetra_SerialComm());
	  }
	else
	  {
	    TEST_FOR_EXCEPTION(true, std::runtime_error, 
			       "Could not find communicator "
			       "for vector space " << tsfSpace);
	  }
      }
    else
      {
	getComm(tsfSpace.getBlock(0), comm);
	int myRank = comm->MyPID();
	int np = comm->NumProc();
	for (int i=1; i<tsfSpace.numBlocks(); i++)
	  {
	    getComm(tsfSpace.getBlock(i), comm);
	    TEST_FOR_EXCEPTION(myRank != comm->MyPID(), std::runtime_error,
			       "inconsistent processor ranks in "
			       "block vector space " << tsfSpace);
	    TEST_FOR_EXCEPTION(np != comm->NumProc(), std::runtime_error,
			       "inconsistent processor counts in "
			       "block vector space " << tsfSpace);
	  }
      }
#else
    comm = rcp(new Epetra_SerialComm());
#endif
  }

}

#undef MPIVectorSpaceBase 





