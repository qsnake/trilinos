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

// We want to solve a 2D finite element problem on the following 2D domain:
//
//          
//   (12)--[21]---(13)--[22]---(14)--[23]--(15)
//     |            |            |           |
//   [17]   <6>   [18]   <7>   [19]   <8>  [20]
//     |            |            |           |
//    (8)--[14]----(9)--[15]---(10)--[16]--(11)
//     |            |            |           |
//   [10]   <3>   [11]   <4>   [12]   <5>  [13]
//     |            |            |           |
//    (4)---[7]----(5)---[8]----(6)---[9]---(7)
//     |            |            |           |
//    [3]   <0>    [4]   <1>    [5]   <2>   [6]
//     |            |            |           |
//    (0)---[0]----(1)---[1]----(2)---[2]---(3)
//
// The notation is:
// *) (X) --> X is a vertex 
// *) <X> --> X is a finite element
// *) [X] --> X is an edge
//
// The grid is subdivided among two processors. We suppose that
// processor 0 owns elements 0, 3, 6 plus 1, 4, 7, while processor 2
// owns elements 2, 5, 8 plus 1, 4, 8. Therefore, elements 1, 4, 7 constitute
// the overlap between the two subdomains
// The total number of unknowns is: 16 + 24 + 9 = 49. These unknowns are 
// distributed as follows, where * indicates processor 0 and @ processor 1.
//
//    (*)---[*]----(*)---[*]----(@)---[@]---(@)
//     |            |            |           |
//    [*]   <*>    [*]   <*>    [@]   <@>   [@]
//     |            |            |           |
//    (*)---[*]----(*)---[*]----(@)---[@]---(@)
//     |            |            |           |
//    [*]   <*>    [*]   <*>    [@]   <@>   [@]
//     |            |            |           |
//    (*)---[*]----(*)---[*]----(@)---[@]---(@)
//     |            |            |           |
//    [*]   <*>    [*]   <*>    [@]   <@>   [@]
//     |            |            |           |
//    (*)---[*]----(*)---[*]----(@)---[@]---(@)
//
// 
// NOTE: The code is very long because I create all the grid structures
// in the code. 
//
// \author Marzio Sala, ETHZ/COLAB
//
// \date Last updated on 29-Nov-05

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
  
  // ============================== //
  // build the grid data structures //
  // ============================== //

  // This is a very simple structure; for each element,
  // I insert the local ID of vertices, then local ID of edges, finally local ID of element.
  // Vertices, edges and elements refer to their local ordering (that is, the
  // first vertex is 0, as the first edge and element). The data refer to the
  // grid with the overlapping region in the middle.
  
  Teuchos::SerialDenseMatrix<int, int> FiniteElements(6, 9);

  if (Comm.getMyImageID() == 0) 
  {
    // vertices
    FiniteElements(0,0) = 0;  FiniteElements(0,1) = 1;  FiniteElements(0,2) = 5;  FiniteElements(0,3) = 4;
    FiniteElements(1,0) = 1;  FiniteElements(1,1) = 2;  FiniteElements(1,2) = 6;  FiniteElements(1,3) = 5;
    FiniteElements(2,0) = 4;  FiniteElements(2,1) = 5;  FiniteElements(2,2) = 9;  FiniteElements(2,3) = 8;
    FiniteElements(3,0) = 5;  FiniteElements(3,1) = 6;  FiniteElements(3,2) = 10; FiniteElements(3,3) = 9;
    FiniteElements(4,0) = 8;  FiniteElements(4,1) = 9;  FiniteElements(4,2) = 13; FiniteElements(4,3) = 12;
    FiniteElements(5,0) = 9;  FiniteElements(5,1) = 10; FiniteElements(5,2) = 14; FiniteElements(5,3) = 13;
    // edges
    FiniteElements(0,4) = 0;  FiniteElements(0,5) = 4;  FiniteElements(0,6) = 7;  FiniteElements(0,7) = 3;
    FiniteElements(1,4) = 1;  FiniteElements(1,5) = 5;  FiniteElements(1,6) = 8;  FiniteElements(1,7) = 4;
    FiniteElements(2,4) = 7;  FiniteElements(2,5) = 11; FiniteElements(2,6) = 14; FiniteElements(2,7) = 10;
    FiniteElements(3,4) = 8;  FiniteElements(3,5) = 12; FiniteElements(3,6) = 15; FiniteElements(3,7) = 11;
    FiniteElements(4,4) = 14; FiniteElements(4,5) = 18; FiniteElements(4,6) = 21; FiniteElements(4,7) = 17;
    FiniteElements(5,4) = 15; FiniteElements(5,5) = 19; FiniteElements(5,6) = 22; FiniteElements(5,7) = 18;
    // elements
    FiniteElements(0,8) = 0;
    FiniteElements(1,8) = 1;
    FiniteElements(2,8) = 3;
    FiniteElements(3,8) = 4;
    FiniteElements(4,8) = 6;
    FiniteElements(5,8) = 7;
  }
  else{
    // vertices
    FiniteElements(0,0) = 1;  FiniteElements(0,1) = 2;  FiniteElements(0,2) = 6;  FiniteElements(0,3) = 5;
    FiniteElements(1,0) = 2;  FiniteElements(1,1) = 3;  FiniteElements(1,2) = 7;  FiniteElements(1,3) = 6;
    FiniteElements(2,0) = 5;  FiniteElements(2,1) = 6;  FiniteElements(2,2) = 10; FiniteElements(2,3) = 9;
    FiniteElements(3,0) = 6;  FiniteElements(3,1) = 7;  FiniteElements(3,2) = 11; FiniteElements(3,3) = 10;
    FiniteElements(4,0) = 9;  FiniteElements(4,1) = 10; FiniteElements(4,2) = 14; FiniteElements(4,3) = 13;
    FiniteElements(5,0) = 10; FiniteElements(5,1) = 11; FiniteElements(5,2) = 15; FiniteElements(5,3) = 14;
    // edges
    FiniteElements(0,4) = 1;  FiniteElements(0,5) = 5;  FiniteElements(0,6) = 8;  FiniteElements(0,7) = 4;
    FiniteElements(1,4) = 2;  FiniteElements(1,5) = 6;  FiniteElements(1,6) = 9;  FiniteElements(1,7) = 5;
    FiniteElements(2,4) = 8;  FiniteElements(2,5) = 12; FiniteElements(2,6) = 15; FiniteElements(2,7) = 11;
    FiniteElements(3,4) = 9;  FiniteElements(3,5) = 13; FiniteElements(3,6) = 16; FiniteElements(3,7) = 12;
    FiniteElements(4,4) = 15; FiniteElements(4,5) = 19; FiniteElements(4,6) = 22; FiniteElements(4,7) = 18;
    FiniteElements(5,4) = 16; FiniteElements(5,5) = 20; FiniteElements(5,6) = 23; FiniteElements(5,7) = 19;
    // elements
    FiniteElements(0,8) = 1;
    FiniteElements(1,8) = 2;
    FiniteElements(2,8) = 4;
    FiniteElements(3,8) = 5;
    FiniteElements(4,8) = 7;
    FiniteElements(5,8) = 8;
  }

  // =================================================== //
  // Creation of spaces related to the grid distribution //
  //                                                     //
  // All variables                                       //
  // have the prefix `Ext' to distinguish between the    //
  // spaces used in the matrix assembly routines. `Ext'  //
  // means the vertices, edges and elements that refer   //
  // to the local part of the grid, while without `Ext'  //
  // we refer to vertices, edges and elements that are   //
  // assigned to the calling image.                      //
  // `Ext' is `Extended' in the sense that we add every- //
  // thing that is included in the overlapping region.   //
  // =================================================== //

  OrdinalType NumMyExtVertices;
  OrdinalType NumGlobalExtVertices = 16;
  vector<OrdinalType> MyGlobalExtVertices;
  
  if (Comm.getMyImageID() == 0)
  {
    NumMyExtVertices = OrdinalOne * 12;
    MyGlobalExtVertices.resize(NumMyExtVertices);
    MyGlobalExtVertices[0] = 0;
    MyGlobalExtVertices[1] = 1;
    MyGlobalExtVertices[2] = 2;
    MyGlobalExtVertices[3] = 4;
    MyGlobalExtVertices[4] = 5;
    MyGlobalExtVertices[5] = 6;
    MyGlobalExtVertices[6] = 8;
    MyGlobalExtVertices[7] = 9;
    MyGlobalExtVertices[8] = 10;
    MyGlobalExtVertices[9] = 12;
    MyGlobalExtVertices[10] = 13;
    MyGlobalExtVertices[11] = 14;
  }
  else
  {
    NumMyExtVertices = OrdinalOne * 12;
    MyGlobalExtVertices.resize(NumMyExtVertices);
    MyGlobalExtVertices[0] = 1;
    MyGlobalExtVertices[1] = 2;
    MyGlobalExtVertices[2] = 3;
    MyGlobalExtVertices[3] = 5;
    MyGlobalExtVertices[4] = 6;
    MyGlobalExtVertices[5] = 7;
    MyGlobalExtVertices[6] = 9;
    MyGlobalExtVertices[7] = 10;
    MyGlobalExtVertices[8] = 11;
    MyGlobalExtVertices[9] = 13;
    MyGlobalExtVertices[10] = 14;
    MyGlobalExtVertices[11] = 15;
  }

  Tpetra::ElementSpace<OrdinalType> ExtVertexSpace(-OrdinalOne, NumMyExtVertices, MyGlobalExtVertices, indexBase, platformE);
  ExtVertexSpace.setLabel("ExtVertesSpace");
  Tpetra::VectorSpace<OrdinalType, ScalarType> VectorExtVertexSpace(ExtVertexSpace, platformV);
  VectorExtVertexSpace.setLabel("ExtVectorVertesSpace");

  // Distribution of edges
  
  OrdinalType NumMyExtEdges;
  OrdinalType NumGlobalEdges = 24;
  vector<OrdinalType> MyGlobalExtEdges;
  
  if (Comm.getMyImageID() == 0)
  {
    NumMyExtEdges = OrdinalOne * 17;
    MyGlobalExtEdges.resize(NumMyExtEdges);
    MyGlobalExtEdges[0] = 0;
    MyGlobalExtEdges[1] = 1;
    MyGlobalExtEdges[2] = 3;
    MyGlobalExtEdges[3] = 4;
    MyGlobalExtEdges[4] = 5;
    MyGlobalExtEdges[5] = 7;
    MyGlobalExtEdges[6] = 8;
    MyGlobalExtEdges[7] = 10;
    MyGlobalExtEdges[8] = 11;
    MyGlobalExtEdges[9] = 12;
    MyGlobalExtEdges[10] = 14;
    MyGlobalExtEdges[11] = 15;
    MyGlobalExtEdges[12] = 17;
    MyGlobalExtEdges[13] = 18;
    MyGlobalExtEdges[14] = 19;
    MyGlobalExtEdges[15] = 21;
    MyGlobalExtEdges[16] = 22;
  }
  else
  {
    NumMyExtEdges = OrdinalOne * 17;
    MyGlobalExtEdges.resize(NumMyExtEdges);
    MyGlobalExtEdges[0] = 1;
    MyGlobalExtEdges[1] = 2;
    MyGlobalExtEdges[2] = 4;
    MyGlobalExtEdges[3] = 5;
    MyGlobalExtEdges[4] = 6;
    MyGlobalExtEdges[5] = 8;
    MyGlobalExtEdges[6] = 9;
    MyGlobalExtEdges[7] = 11;
    MyGlobalExtEdges[8] = 12;
    MyGlobalExtEdges[9] = 13;
    MyGlobalExtEdges[10] = 15;
    MyGlobalExtEdges[11] = 16;
    MyGlobalExtEdges[12] = 18;
    MyGlobalExtEdges[13] = 19;
    MyGlobalExtEdges[14] = 20;
    MyGlobalExtEdges[15] = 22;
    MyGlobalExtEdges[16] = 23;
  }

  Tpetra::ElementSpace<OrdinalType> ExtEdgeSpace(-OrdinalOne, NumMyExtEdges, MyGlobalExtEdges, indexBase, platformE);
  ExtEdgeSpace.setLabel("ExtEdgeSpace");
  Tpetra::VectorSpace<OrdinalType, ScalarType> VectorExtEdgeSpace(ExtEdgeSpace, platformV);
  VectorExtEdgeSpace.setLabel("VectorExtEdgeSpace");

  // Distribution of elements
  
  OrdinalType NumMyExtFiniteElements;
  OrdinalType NumGlobalFiniteElements = 9;
  vector<OrdinalType> MyGlobalExtFiniteElements;
  
  if (Comm.getMyImageID() == 0)
  {
    NumMyExtFiniteElements = OrdinalOne * 6;
    MyGlobalExtFiniteElements.resize(NumMyExtFiniteElements);
    MyGlobalExtFiniteElements[0] = 0;
    MyGlobalExtFiniteElements[1] = 1;
    MyGlobalExtFiniteElements[2] = 3;
    MyGlobalExtFiniteElements[3] = 4;
    MyGlobalExtFiniteElements[4] = 6;
    MyGlobalExtFiniteElements[5] = 7;
  }
  else
  {
    NumMyExtFiniteElements = OrdinalOne * 6;
    MyGlobalExtFiniteElements.resize(NumMyExtFiniteElements);
    MyGlobalExtFiniteElements[0] = 1;
    MyGlobalExtFiniteElements[1] = 2;
    MyGlobalExtFiniteElements[2] = 4;
    MyGlobalExtFiniteElements[3] = 5;
    MyGlobalExtFiniteElements[4] = 7;
    MyGlobalExtFiniteElements[5] = 8;
  }

  Tpetra::ElementSpace<OrdinalType> ExtFiniteElementSpace(-OrdinalOne, NumMyExtFiniteElements, MyGlobalExtFiniteElements, indexBase, platformE);
  ExtFiniteElementSpace.setLabel("ExtFiniteElementSpace");
  Tpetra::VectorSpace<OrdinalType, ScalarType> VectorExtFiniteElementSpace(ExtFiniteElementSpace, platformV);
  VectorExtFiniteElementSpace.setLabel("VectorExtFiniteElementSpace");

  // ================================== //
  // Creation of non-overlapping spaces //
  // ================================== //
  
  // Distribution of vertices
  
  OrdinalType NumMyVertices;
  OrdinalType NumGlobalVertices = 16;
  vector<OrdinalType> MyGlobalVertices;
  
  if (Comm.getMyImageID() == 0)
  {
    NumMyVertices = OrdinalOne * 8;
    MyGlobalVertices.resize(NumMyVertices);
    MyGlobalVertices[0] = 0;
    MyGlobalVertices[1] = 1;
    MyGlobalVertices[2] = 4;
    MyGlobalVertices[3] = 5;
    MyGlobalVertices[4] = 8;
    MyGlobalVertices[5] = 9;
    MyGlobalVertices[6] = 12;
    MyGlobalVertices[7] = 13;
  }
  else
  {
    NumMyVertices = OrdinalOne * 8;
    MyGlobalVertices.resize(NumMyVertices);
    MyGlobalVertices[0] = 2;
    MyGlobalVertices[1] = 3;
    MyGlobalVertices[2] = 6;
    MyGlobalVertices[3] = 7;
    MyGlobalVertices[4] = 10;
    MyGlobalVertices[5] = 11;
    MyGlobalVertices[6] = 14;
    MyGlobalVertices[7] = 15;
  }

  Tpetra::ElementSpace<OrdinalType> VertexSpace(-OrdinalOne, NumMyVertices, MyGlobalVertices, indexBase, platformE);
  VertexSpace.setLabel("VertexSpace");
  Tpetra::VectorSpace<OrdinalType, ScalarType> VectorVertexSpace(VertexSpace, platformV);
  VectorVertexSpace.setLabel("VectorVertexSpace");

  // Distribution of edges
  
  OrdinalType NumMyEdges;
  vector<OrdinalType> MyGlobalEdges;
  
  if (Comm.getMyImageID() == 0)
  {
    NumMyEdges = OrdinalOne * 14;
    MyGlobalEdges.resize(NumMyEdges);
    MyGlobalEdges[0] = 0;
    MyGlobalEdges[1] = 1;
    MyGlobalEdges[2] = 3;
    MyGlobalEdges[3] = 4;
    MyGlobalEdges[4] = 7;
    MyGlobalEdges[5] = 8;
    MyGlobalEdges[6] = 10;
    MyGlobalEdges[7] = 11;
    MyGlobalEdges[8] = 14;
    MyGlobalEdges[9] = 15;
    MyGlobalEdges[10] = 17;
    MyGlobalEdges[11] = 18;
    MyGlobalEdges[12] = 21;
    MyGlobalEdges[13] = 22;
  }
  else
  {
    NumMyEdges = OrdinalOne * 10;
    MyGlobalEdges.resize(NumMyEdges);
    MyGlobalEdges[0] = 2;
    MyGlobalEdges[1] = 5;
    MyGlobalEdges[2] = 6;
    MyGlobalEdges[3] = 9;
    MyGlobalEdges[4] = 12;
    MyGlobalEdges[5] = 13;
    MyGlobalEdges[6] = 16;
    MyGlobalEdges[7] = 19;
    MyGlobalEdges[8] = 20;
    MyGlobalEdges[9] = 23;
  }

  Tpetra::ElementSpace<OrdinalType> EdgeSpace(-OrdinalOne, NumMyEdges, MyGlobalEdges, indexBase, platformE);
  EdgeSpace.setLabel("EdgeSpace");
  Tpetra::VectorSpace<OrdinalType, ScalarType> VectorEdgeSpace(EdgeSpace, platformV);
  VectorEdgeSpace.setLabel("VertexEdgeSpace");

  // Distribution of elements
  
  OrdinalType NumMyFiniteElements;
  vector<OrdinalType> MyGlobalFiniteElements;
  
  if (Comm.getMyImageID() == 0)
  {
    NumMyFiniteElements = OrdinalOne * 6;
    MyGlobalFiniteElements.resize(NumMyFiniteElements);
    MyGlobalFiniteElements[0] = 0;
    MyGlobalFiniteElements[1] = 1;
    MyGlobalFiniteElements[2] = 3;
    MyGlobalFiniteElements[3] = 4;
    MyGlobalFiniteElements[4] = 6;
    MyGlobalFiniteElements[5] = 7;
  }
  else
  {
    NumMyFiniteElements = OrdinalOne * 3;
    MyGlobalFiniteElements.resize(NumMyFiniteElements);
    MyGlobalFiniteElements[0] = 2;
    MyGlobalFiniteElements[1] = 5;
    MyGlobalFiniteElements[2] = 8;
  }

  Tpetra::ElementSpace<OrdinalType> FiniteElementSpace(-OrdinalOne, NumMyFiniteElements, MyGlobalFiniteElements, indexBase, platformE);
  FiniteElementSpace.setLabel("FiniteElementSpace");
  Tpetra::VectorSpace<OrdinalType, ScalarType> VectorFiniteElementSpace(FiniteElementSpace, platformV);
  VectorFiniteElementSpace.setLabel("VectorFiniteElementSpace");

  // ==================================================================== //
  // At this point we have the distribution of vertices, edges and        //
  // elements across the processors. We can compose the "final" map,      //
  // which contains the distribution of unknowns for the final matrix.    //
  // We suppose to order vertices first, followed by edges, then          //
  // elements. For this purpose we build another couple of spaces         //
  // (RowSpace and VectorRowSpace), which will contain this distribution. //
  // ==================================================================== //

  int NumMyRows = NumMyVertices + NumMyEdges + NumMyFiniteElements;
  vector<OrdinalType> MyGlobalRows(NumMyRows);

  OrdinalType count = OrdinalZero;

  // this are the global ID's of vertices
  for (OrdinalType i = OrdinalZero ; i < NumMyVertices ; ++i)
  {
    MyGlobalRows[count++] = MyGlobalVertices[i];
  }

  // this are the global ID's of edges
  for (OrdinalType i = OrdinalZero ; i < NumMyEdges ; ++i)
  {
    MyGlobalRows[count++] = MyGlobalEdges[i] + NumGlobalVertices;
  }

  // this are the global ID's of finite elements
  for (OrdinalType i = OrdinalZero ; i < NumMyFiniteElements ; ++i)
  {
    MyGlobalRows[count++] = MyGlobalFiniteElements[i] + NumGlobalVertices + NumGlobalEdges;
  }

  assert (count == NumMyRows);

  Tpetra::ElementSpace<OrdinalType> RowSpace(-OrdinalOne, NumMyRows, MyGlobalRows, indexBase, platformE);
  Tpetra::VectorSpace<OrdinalType, ScalarType> VectorRowSpace(RowSpace, platformV);
  VectorRowSpace.setLabel("VectorRowSpace");
  
  // ================ //
  // Setup the matrix //
  // ================ //
  
  OrdinalType NumUnknownsPerFiniteElement = 4 + 4 + 1;
  Teuchos::SerialDenseMatrix<OrdinalType, ScalarType> localMatrix(NumUnknownsPerFiniteElement, NumUnknownsPerFiniteElement);
  localMatrix.putScalar(ScalarZero);
  // diagonal now 
  for (OrdinalType i = OrdinalZero ; i < NumUnknownsPerFiniteElement ; ++i)
    localMatrix(i, i) = ScalarOne;

  // Allocate the matrix, based on VectorRowSpace
  
  Tpetra::CisMatrix<OrdinalType,ScalarType> matrix(VectorRowSpace);

  // add zero diagonal

  for (OrdinalType LID = OrdinalZero ; LID < NumMyRows ; ++LID)
  {
    OrdinalType GID = MyGlobalRows[LID];
    matrix.submitEntry(Tpetra::Insert, GID, ScalarZero, GID);
  }

  // Loop over all the (locally owned) elements
  
  for (OrdinalType FEID = OrdinalZero ; FEID < NumMyExtFiniteElements ; ++FEID)
  {
    vector<OrdinalType> LIDs(NumUnknownsPerFiniteElement);
    vector<OrdinalType> GIDs(NumUnknownsPerFiniteElement);

    // load the global IDs for vertices, edges and elements

    for (OrdinalType i = OrdinalZero ; i < OrdinalOne * 4 ; ++i)
    {
      GIDs[i] = FiniteElements(FEID, i);
      LIDs[i] = ExtVertexSpace.getLID(GIDs[i]);
    }

    for (OrdinalType i = OrdinalZero ; i < OrdinalOne * 4 ; ++i)
    {
      GIDs[i + 4] = FiniteElements(FEID, i + 4) + NumGlobalVertices;
      LIDs[i + 4] = ExtEdgeSpace.getLID(FiniteElements(FEID, i + 4)) + NumGlobalVertices;
    }

    GIDs[8] = FiniteElements(FEID, 8) + NumGlobalVertices + NumGlobalEdges;
    LIDs[8] = ExtFiniteElementSpace.getLID(FiniteElements(FEID, 8)) + NumGlobalVertices + NumGlobalEdges;

    // get the coordinates (NOT DONE HERE)

    // build the local matrix, in this case A_loc (NOT DONE HERE)

    // submit entries of localMatrix into the matrix

    for (OrdinalType row = OrdinalZero ; row < NumUnknownsPerFiniteElement ; ++row)
    {
      // We skip non-locally owned vertices
      if (RowSpace.isMyGID(GIDs[row]))
      {
        for (OrdinalType col = OrdinalZero ; col < NumUnknownsPerFiniteElement ; ++col)
        {
          matrix.submitEntry(Tpetra::Add, GIDs[row], localMatrix(row, col), GIDs[col]);
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
