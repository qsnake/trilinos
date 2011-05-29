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
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_ScalarTraits.hpp"

typedef int OrdinalType;
typedef double ScalarType;

// We want to solve the 1D finite element problem
//
//    - u'' = f   x \in \Omega
//      u(0) = u_0
//      u(1) = u_1
//
// on the following domain defined by 5 elements and 6 nodes:
//
//       e0    e1    e2    e3    e4  
//     o --- o --- o --- o --- o --- o 
//     ^                             ^
//     n0    n1    n2    n3    n4    n5
//
//     |<--------->|     |<--------->|
//        proc 0             proc 1
//
// The grid is subdivided as follows among the two processors:
// (where (*) defines a locally owned node, and [*] a ghost node)
//
// proc 0:
//
//     e0     e1     e2
// (0) -- (1) -- (2) -- [3] -- (X) -- (X)
//
// proc 1:
//
//                   e2     e0     e1
// (X) -- (X) -- [3] -- (0) -- (1) -- (2)
//
// \author Marzio Sala, ETHZ/COLAB
//
// \date Last updated on 28-Nov-05

// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[]) 
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Tpetra::MpiComm<OrdinalType, ScalarType> Comm(MPI_COMM_WORLD);
#else
  Tpetra::SerialComm<OrdinalType, ScalarType> Comm;
#endif

  if (Comm.getNumImages() != 2)
  {
    cout << "This example must be ran with two processors" << endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_SUCCESS);
  }

  // Get zero and one for the OrdinalType
  
  OrdinalType const OrdinalZero = Teuchos::ScalarTraits<OrdinalType>::zero();
  OrdinalType const OrdinalOne  = Teuchos::ScalarTraits<OrdinalType>::one();

  // Get zero and one for the ScalarType
  
  ScalarType const ScalarOne  = Teuchos::ScalarTraits<ScalarType>::one();
  ScalarType const ScalarZero = Teuchos::ScalarTraits<ScalarType>::zero();

  OrdinalType indexBase = OrdinalZero;

  // Creation of a platform
  
#ifdef HAVE_MPI
  const Tpetra::MpiPlatform <OrdinalType, OrdinalType> platformE(MPI_COMM_WORLD);
  const Tpetra::MpiPlatform <OrdinalType, ScalarType> platformV(MPI_COMM_WORLD);
#else
  const Tpetra::SerialPlatform <OrdinalType, OrdinalType> platformE;
  const Tpetra::SerialPlatform <OrdinalType, ScalarType> platformV;
#endif

  // These are the `elements' assigned to this image
  // NOTE: This element is NOT a finite element, instead it is
  //       a component of the ElementSpace (that is, the distribution of
  //       vertices). FiniteElements values are defined later on.
  // 
  // NOTATION:
  // - Element       -> an entity of the space
  // - Vertex        -> finite element vertex
  // - FiniteElement -> finite element, composed by two Vertices
  
  const OrdinalType NumMyVertices = OrdinalOne * 3;
  vector<OrdinalType> MyGlobalVertices(NumMyVertices);

  if (Comm.getMyImageID() == 0)
  {
    for (OrdinalType i = OrdinalZero ; i < NumMyVertices ; ++i)
      MyGlobalVertices[OrdinalZero + i] = OrdinalZero + i;
  }
  else
  {
    for (OrdinalType i = OrdinalZero ; i < NumMyVertices ; ++i)
      MyGlobalVertices[OrdinalZero + i] = OrdinalZero + i + 3;
  }

  Tpetra::ElementSpace<OrdinalType> VertexSpace(-OrdinalOne, NumMyVertices, MyGlobalVertices, indexBase, platformE);
  Tpetra::VectorSpace<OrdinalType, ScalarType> VectorVertexSpace(VertexSpace, platformV);

  const OrdinalType NumMyPaddedVertices = OrdinalOne * 4;
  vector<OrdinalType> MyGlobalPaddedVertices(NumMyPaddedVertices);
  // Vector BoundaryVertices contains a flag that specifies
  // whether the node is a Dirichlet node or not
  vector<bool> BoundaryVertices(NumMyVertices);

  // Setting the elements of the connectivity
  if (Comm.getMyImageID() == 0)
  {
    MyGlobalPaddedVertices[0] = OrdinalZero;
    MyGlobalPaddedVertices[1] = OrdinalOne;
    MyGlobalPaddedVertices[2] = OrdinalOne * 2;
    MyGlobalPaddedVertices[3] = OrdinalOne * 3; // ghost node

    BoundaryVertices[0] = true;
    BoundaryVertices[1] = false;
    BoundaryVertices[2] = false;
  }
  else
  {
    MyGlobalPaddedVertices[0] = OrdinalOne * 3;
    MyGlobalPaddedVertices[1] = OrdinalOne * 4;
    MyGlobalPaddedVertices[2] = OrdinalOne * 5;
    MyGlobalPaddedVertices[3] = OrdinalOne * 2; // ghost node

    BoundaryVertices[0] = false;
    BoundaryVertices[1] = false;
    BoundaryVertices[2] = true;
  }

  // This space is needed to import the values of coordinates corresponding
  // to ghost elements. 
  Tpetra::ElementSpace<OrdinalType> PaddedVertexSpace(-OrdinalOne, NumMyPaddedVertices, MyGlobalPaddedVertices, indexBase, platformE);
  Tpetra::VectorSpace<OrdinalType, ScalarType> VectorPaddedVertexSpace(PaddedVertexSpace, platformV);

  // This is the connectivity, in local numbering

  const OrdinalType NumVerticesPerFiniteElement = 2;
  const OrdinalType NumDimensions = 2;

  const OrdinalType NumMyElements = OrdinalOne * 3;
  Teuchos::SerialDenseMatrix<OrdinalType, OrdinalType> Connectivity(NumMyElements, NumVerticesPerFiniteElement);

  // this contains the coordinates
  Tpetra::Vector<OrdinalType, ScalarType> Coord(VectorPaddedVertexSpace);

  if (Comm.getMyImageID() == 0)
  {
    Connectivity(0,0) = 0; Connectivity(0,1) = 1;
    Connectivity(1,0) = 1; Connectivity(1,1) = 2;
    Connectivity(2,0) = 2; Connectivity(2,1) = 3;
    Coord[0] = ScalarZero;
    Coord[1] = ScalarOne;
    Coord[2] = ScalarOne * 2;
    Coord[3] = ScalarOne * 3; // ghost node is last
  }
  else
  {
    Connectivity(0,0) = 0; Connectivity(0,1) = 1;
    Connectivity(1,0) = 1; Connectivity(1,1) = 2;
    Connectivity(2,0) = 3; Connectivity(2,1) = 0;
    Coord[0] = ScalarOne * 3;
    Coord[1] = ScalarOne * 4;
    Coord[2] = ScalarOne * 5;
    Coord[3] = ScalarOne * 2; // ghost node is last
  }

  // Setup the matrix
  
  Teuchos::SerialDenseMatrix<OrdinalType, ScalarType> localMatrix(NumVerticesPerFiniteElement, NumVerticesPerFiniteElement);

  Tpetra::CisMatrix<OrdinalType,ScalarType> matrix(VectorVertexSpace);

  // Loop over all the (locally owned) elements
  
  for (OrdinalType FEID = OrdinalZero ; FEID < NumMyElements ; ++FEID)
  {
    vector<OrdinalType> LIDs(NumVerticesPerFiniteElement), GIDs(NumVerticesPerFiniteElement);

    // Get the local and global ID of this element's vertices
    for (OrdinalType i = OrdinalZero ; i < NumVerticesPerFiniteElement ; ++i)
    {
      LIDs[i] = Connectivity(FEID, i);
      GIDs[i] = PaddedVertexSpace.getGID(LIDs[i]);
    }

    // get the coordinates of all nodes in the element
    vector<ScalarType> x(NumVerticesPerFiniteElement);

    for (OrdinalType i = 0 ; i < NumDimensions ; ++i)
    {
      x[i] = Coord[LIDs[i]];
    }
    
    // build the local matrix, in this case A_loc;
    // this is specific to this 1D Laplace, and does not use
    // coordinates. 
    localMatrix(0, 0) = ScalarOne;
    localMatrix(1, 0) = -ScalarOne;
    localMatrix(0, 1) = -ScalarOne;
    localMatrix(1, 1) = ScalarOne;

    // submit entries of localMatrix into the matrix

    for (OrdinalType LRID = OrdinalZero ; LRID < NumVerticesPerFiniteElement ; ++LRID)
    {
      OrdinalType LocalRow = LIDs[LRID];

      // We skip non-locally owned vertices
      if (LocalRow < NumMyVertices)
      {
        // If the node is a boundary node, then put 1 on the diagonal
        if (BoundaryVertices[LocalRow])
        {
          matrix.submitEntry(Tpetra::Insert, GIDs[LRID], ScalarOne, GIDs[LRID]);
        }
        else
        {
          // otherwise add the entire row
          for (OrdinalType LCID = OrdinalZero ; LCID < NumVerticesPerFiniteElement ; ++LCID)
          {
            matrix.submitEntry(Tpetra::Add, GIDs[LRID], localMatrix(LRID,LCID), GIDs[LCID]);
          }
        }
      }
    }
  }

  matrix.fillComplete();

  cout << matrix;
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return(EXIT_SUCCESS);
}
