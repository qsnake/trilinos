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

#include "TSFEpetraVectorType.hpp"
#include "TSFEpetraVectorSpace.hpp"
#include "TSFEpetraGhostImporter.hpp"
#include "TSFEpetraMatrixFactory.hpp"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Comm.h"
#include "SundanceOut.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Teuchos_RefCountPtr.hpp"
#include "TSFEpetraMatrix.hpp"

using namespace TSFExtended;
using namespace Teuchos;

EpetraVectorType::EpetraVectorType()
{;}


RCP<const Thyra::VectorSpaceBase<double> > 
EpetraVectorType::createSpace(int /*dimension*/,
  int nLocal,
  const int* localIndices,
  const MPIComm& comm) const
{
#ifdef HAVE_MPI
  Epetra_MpiComm epComm(comm.getComm());
#else
  Epetra_SerialComm epComm;
#endif

  TEST_FOR_EXCEPTION(nLocal < 0, std::runtime_error, "negative vector size n=" << nLocal);

	RCP<Epetra_Map> map = rcp(new Epetra_Map(-1, nLocal,
      (int*) localIndices,
      0, epComm));

	return rcp(new EpetraVectorSpace(map));
}

RCP<GhostImporter<double> > 
EpetraVectorType::createGhostImporter(const VectorSpace<double>& space,
                                      int nGhost,
                                      const int* ghostIndices) const
{
  const EpetraVectorSpace* p 
    = dynamic_cast<const EpetraVectorSpace*>(space.ptr().get());

  TEST_FOR_EXCEPTION(p==0, std::runtime_error,
                     "non-epetra vector space [" << space.description() << "] given as "
                     "argument to EpetraVectorType::createGhostImporter()");

  return rcp(new EpetraGhostImporter(p->epetraMap(), nGhost, ghostIndices));
  
}

RCP<MatrixFactory<double> >
EpetraVectorType::createMatrixFactory(const VectorSpace<double>& domain,
                                      const VectorSpace<double>& range) const
{
  RCP<const EpetraVectorSpace> pd 
    = rcp_dynamic_cast<const EpetraVectorSpace>(domain.ptr());

  RCP<const EpetraVectorSpace> pr 
    = rcp_dynamic_cast<const EpetraVectorSpace>(range.ptr());


  TEST_FOR_EXCEPTION(pd.get()==0, std::runtime_error, 
                     "incompatible domain space given to "
                     "EpetraVectorType::createMatrix()");

  TEST_FOR_EXCEPTION(pr.get()==0, std::runtime_error, 
                     "incompatible range space given to "
                     "EpetraVectorType::createMatrix()");

  //  RCP<SingleScalarTypeOp<double> > A = rcp(new EpetraMatrix(pd, pr));

  return rcp(new EpetraMatrixFactory(pd, pr));
}






