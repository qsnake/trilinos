
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

// -----------------
// The files in this directory solve the Chandrasekhar H-equation using NOX. 
// The specific examples that use this file:  ex3.cpp, ex4.cpp, and ex5.cpp.
//
// The Chandrasekhar H-equation is defined as:
//   F(H)(u) = H(u) - [ 1 - c/2 * Int_0^1 { u*H(v)dv / (u+v) } ]^{-1}
// The discretized H-equation is:
//   F(x)_i = x_i - [ 1 - c/(2N) * Sum_{j=1}^N { u_i*x_j / (u_i+u_j) } ]^{-1}
// where u_i = (i - 0.5)/N for 1 <= i <= N.
// 
// The H-equation has two solutions for c in (0,1).  It has one solution for c=1,
//   which is also the turning point for the solution graph.
// 
// This file contains the problem class and the definitions of 
// functions and variables associated with it.
//
 
#ifndef HEQ_H
#define HEQ_H

// Below are all the necessary include files.

// Used for the overlap map
#include "Epetra_Import.h"
#include "Epetra_Export.h"
// Used for output
#include <iostream>
// check for parallel computation
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
// Epetra definition files
#include "Epetra_ConfigDefs.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
// Nox definition files
#include "NOX.H"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_LinearSystem_AztecOO.H"
#include "NOX_Epetra_Group.H"
// Loca definition files
#include "LOCA_Parameter_Vector.H"
#include "LOCA_Epetra_Interface_TimeDependent.H"

// Set up the Problem class that defines the H-equation.  All of the specific
// code related to solving the H-equation should appear in Heq.H and Heq.C.

class HeqProblem {

public:

// The constructor. Requires the number of grid points N_, 
// the Epetra_Communicator, and the output file name.

  HeqProblem(const int N_, const Epetra_Comm * Comm_, ofstream& file);

// destructor
  ~HeqProblem();

// Compute the value of F(x)
  void ComputeF(const Epetra_Vector & x, Epetra_Vector & f, const double c); 
    
// Compute the Jacobian matrix of F(x) for a given x
  void ComputeJacobian(const Epetra_Vector & x, const double c, double alpha); 

// Creates the map for solution vectors and the Jacobian matrix  
  Epetra_CrsMatrix * CreateH();

// Returns a pointer to the internally stored Jacobian matrix
  Epetra_CrsMatrix * GetMatrix();

// Calculate and print the continuation parameter value and the 1-norm of the 
// solution vector
  void printSolution(const Epetra_Vector &x, const double conParam);

// Returns the Map associated with solution vectors / the Jacobian matrix
  Epetra_Map GetMap();


// Definitions of variables used throughout the Problem class
private:

  int N;                           // number of grid points
  double * hold;                   // pointer to temporary variable (used for calculations)
  const Epetra_Comm * Comm;        // the epetra communicator
  Epetra_CrsMatrix * Matrix;       // the Jacobian matrix
  Epetra_Map Map, OverlapMap;      // Map = map for the solution vector;
                                   // OverlapMap = map used for F(x), Jacobian calculations
                                   //    (due to parallel processing)
  ofstream* outputFilePtr;         // pointer to output file
  bool FirstTime;                  // boolean variable used to determine whether the Jacobian
                                   //   matrix has been accessed previously

}; /* class HeqProblem */
#endif
