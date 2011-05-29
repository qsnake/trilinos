
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
// This file contains the problem interface class and the definitions of 
// functions and variables associated with it.
 
#ifndef PROBINTF_H
#define PROBINTF_H


#include "Heq.H"

// define the Problem Interface class (called SimpleProblemInterface)

class SimpleProblemInterface : 
                               public LOCA::Epetra::Interface::Required,
//  the commented commands below are additional options for future use
//                               public LOCA::Epetra::Interface::TimeDependent,
//                               public LOCA::Epetra::Interface::TimeDependentMatrixFree,
                               public NOX::Epetra::Interface::Jacobian
{

public:
 
  // Constructor
  SimpleProblemInterface(HeqProblem * Problem_, double c_);

  // Destructor
  ~SimpleProblemInterface(); 

  // Compute F(x)
  bool computeF(const Epetra_Vector & x, Epetra_Vector & f,
                NOX::Epetra::Interface::Required::FillType F );
  
  // Compute the Jacobian of F(x)
  bool computeJacobian(const Epetra_Vector & x, Epetra_Operator & Jac);

  bool computePrecMatrix(const Epetra_Vector & x, Epetra_RowMatrix & M);
  
  bool computePreconditioner(const Epetra_Vector & x, Epetra_Operator & O); 

  // Set the parameter vector
  void setParameters(const LOCA::ParameterVector& params);

  // Print out the solution in a specific format
  void printSolution(const Epetra_Vector &x, const double conParam);


  // variable definitions

  HeqProblem * Problem; 
  double c;
  
};
#endif
