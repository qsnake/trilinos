// @HEADER
// ***********************************************************************
// 
//                      Didasko Tutorial Package
//                 Copyright (2005) Sandia Corporation
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
//
// Questions about Didasko? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
// 
// ***********************************************************************
// @HEADER

#include "Tpetra_ConfigDefs.hpp"
#ifdef HAVE_MPI
#include "Tpetra_MpiPlatform.hpp"
#include "Tpetra_MpiComm.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_SerialComm.hpp"
#endif
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_CisMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Teuchos_ScalarTraits.hpp"

// \author Marzio Sala, ETHZ/COLAB
//
// \date Last updated on 28-Nov-05

typedef int OrdinalType;
typedef double ScalarType;

int main(int argc, char *argv[]) 
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Tpetra::MpiComm<OrdinalType, ScalarType> Comm(MPI_COMM_WORLD);
#else
  Tpetra::SerialComm<OrdinalType, ScalarType> Comm;
#endif

  // Get zero and one for the OrdinalType
  
  OrdinalType const OrdinalZero = Teuchos::ScalarTraits<OrdinalType>::zero();
  OrdinalType const OrdinalOne  = Teuchos::ScalarTraits<OrdinalType>::one();

  // Get zero and one for the ScalarType
  
  ScalarType const ScalarZero = Teuchos::ScalarTraits<ScalarType>::zero();
  ScalarType const ScalarOne  = Teuchos::ScalarTraits<ScalarType>::one();

  // Creates a vector of size `length', then set the elements values.
  
  OrdinalType length    = OrdinalOne * 10;
  OrdinalType indexBase = OrdinalZero;

  // 1) Creation of a platform
  
#ifdef HAVE_MPI
  const Tpetra::MpiPlatform <OrdinalType, OrdinalType> platformE(MPI_COMM_WORLD);
  const Tpetra::MpiPlatform <OrdinalType, ScalarType> platformV(MPI_COMM_WORLD);
#else
  const Tpetra::SerialPlatform <OrdinalType, OrdinalType> platformE;
  const Tpetra::SerialPlatform <OrdinalType, ScalarType> platformV;
#endif

  // 2) We can now create a space:

  Tpetra::ElementSpace<OrdinalType> elementSpace(length, indexBase, platformE);
  Tpetra::VectorSpace<OrdinalType, ScalarType> 
    vectorSpace(elementSpace, platformV);

  // 3) and the vector, which has type int for the OrdinalType
  //    and double for the ScalarType

  Tpetra::Vector<OrdinalType, ScalarType> x(vectorSpace), y(vectorSpace);

  x.setAllToRandom();
  y.setAllToScalar(ScalarZero);

  // 4) We now setup the matrix, which is diagonal
  
  OrdinalType NumMyElements = elementSpace.getNumMyElements();
  vector<OrdinalType> MyGlobalElements = elementSpace.getMyGlobalElements();
  
  Tpetra::CisMatrix<OrdinalType,ScalarType> matrix(vectorSpace);

  for (OrdinalType LID = OrdinalZero ; LID < NumMyElements ; ++LID)
  {
    OrdinalType GID = MyGlobalElements[LID];
    // add the diagonal element of value `GID'
    matrix.submitEntry(Tpetra::Insert, GID, (ScalarType)GID, GID);
  }

  matrix.fillComplete();

  // 5) We multiply the matrix for vector x, the result is in y

  matrix.apply(x, y, false);

  ScalarType Norm2 = y.norm2();

  if (Comm.getMyImageID() == 0)
    cout << "|| A * x ||_2 = " << Norm2 << endl;

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return(EXIT_SUCCESS);
}
