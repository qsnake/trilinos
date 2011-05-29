
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
// This file contains the code to compute the functions associated with 
// the problem class.
 

#include "Heq.H"

// The constructor. Requires the number of grid points N, the Epetra communicator Comm,
//   and the output file pointer "file".  
// Also initializes the variables used by the HeqProblem, creates the Map & OverlapMap,
//   and sets the boolean variable FirstTime to 'true'

  HeqProblem::HeqProblem(const int N_, const Epetra_Comm * Comm_, ofstream& file) :
    N(N_), Comm(Comm_), Map(N_,0,*Comm_), OverlapMap(N_,N_,0,*Comm_), outputFilePtr(&file), FirstTime(true)
  {
    Matrix = CreateH();
  }

// Destructor
  HeqProblem::~HeqProblem() 
  {
	 delete Matrix;
  }

//
// Compute the value of F(x)
//
  void HeqProblem::ComputeF(const Epetra_Vector & x, Epetra_Vector & f, const double c) 
    {
    // Set up variables
    int ilocal, j, iglobal;  // counter variables used for loops:
                             //   ilocal = counter for local elements on this processor;
                             //   iglobal = counter to signify global position across all procs 
    int Myele;                // holds the number of elements on the processor
    double temp;             // temporary variable used for calculations
    hold = new double[N];    // array to hold calculated values used to compute F(x) and 
                             //    the Jacobian matrix
  
    // Set up the importer for the overlap map
    Epetra_Import Importer(OverlapMap, Matrix->Map()); 

    // Define the overlap map
    Epetra_Vector overlap_x(OverlapMap);
 
    // Fill in the overlap map
    overlap_x.Import(x, Importer, Insert);

    // Get the no. of elements on this processor
    Myele = Map.NumMyElements();

    // Compute the function values
    //   Loop over the number of elements on the local processor
    for (ilocal=0; ilocal<Myele; ilocal++)  {
      temp=0.0;
      iglobal=Map.GID(ilocal);    // Get the global ID for the local element

      for (j=0; j<N; j++) {       // Compute the value of Sum_{j=1}^N { u_i*x_j / (u_i+u_j) }    
	 		temp = temp + ((double)iglobal+0.5)*overlap_x[j]/(double)(iglobal+j+1);
      }      
      temp = 1.0 - c*temp/(double)(2*N);

    // Store the value of [1-c/2N * Sum()]^{-1} for use in computing the Jacobian 
      hold[iglobal] = pow(temp,-1.0);
      f[ilocal] = overlap_x[iglobal]-hold[iglobal];  // compute final F(x)_i value

    }  // end of ilocal loop

  }  // end computeF

//
// Compute the Jacobian matrix for a given vector x
//
// The discretized version of the Jacobian is:
//   df_i/dx_j = \delta_ij - c/2N * [ u_i / (u_i+u_j) ] * hold[i]^2
// where hold[i] is calculated by the computeF function == Sum_{j=1}^N { u_i*x_j / (u_i+u_j) }  
//
// The Jacobian is calculated row by row and then the values are inserted into the matrix
//
// The function inputs are the solution vector x, the value of c, and a variable alpha that 
//   may be used in future code to compute the value of the shifted matrix, alpha*J + beta*M,
//   where M is the mass matrix.  To get just the Jacobian J, alpha = 1.0

  void HeqProblem::ComputeJacobian(const Epetra_Vector & x, const double c, double alpha) 
  {    
    // Set up variables
    double *Values = new double[N]; // holds the values computed for each row
    int *Indices = new int[N]; // holds the column index corresponding to each calculated value   
    int ilocal, iglobal, j;     // counter variables used for loops:
                               //   ilocal = counter for local elements on this processor;
                               //   iglobal = counter to signify global position across all procs 
    int Myele;                 // holds the number of elements on the processor
    double temp;               // temporary variable used for calculations

    // Set up the importer for the overlap map
    Epetra_Import Importer(OverlapMap, Matrix->Map()); 

    // Define the overlap map
    Epetra_Vector overlap_x(OverlapMap);

    // Fill in the overlap map
    overlap_x.Import(x, Importer, Insert);

    // Get the no. of elements on this processor
    Myele = Map.NumMyElements();

    // Compute the values for the Jacobian matrix
    //   Loop over the number of elements on the local processor
    for (ilocal=0; ilocal<Myele; ilocal++) {
       iglobal=Map.GID(ilocal);   // Get the global ID for the local element
 
       for (j=0; j<N; j++) {      // compute each Jacobian element
         temp = pow(hold[iglobal],2.0);
         temp = -temp*(c/(double)(2*N))*((double)iglobal+0.5)/(double)(iglobal+j+1);			

         if (iglobal==j) 
         	temp = temp + 1.0;
			Values[j] = alpha*temp;  // save the computed value
			Indices[j] = j;          // assign the column index of the computed value

      }   // end of j loop

    // Insert the values calculated into the matrix by row
    if (FirstTime == true) 
      Matrix->InsertGlobalValues(iglobal, N, &Values[0], &Indices[0]);
    else
      Matrix->ReplaceGlobalValues(iglobal, N, &Values[0], &Indices[0]);

  }  // end of ilocal loop

  // Done with matrix calculations!  
  // Do final steps to deallocate memory and finish the matrix fill.
  if (FirstTime == true)
    Matrix->FillComplete();

  delete [] Indices;
  delete [] Values;

  FirstTime = false;
  }  // end of ComputeJacobian

//
// Creates the map for vectors and the Jacobian matrix 
//
  Epetra_CrsMatrix * HeqProblem::CreateH()
{
  // create Jacobian matrix
  Matrix = new Epetra_CrsMatrix(Copy,Map,N);

  return Matrix;
}  

//        
// Returns a pointer to the Jacobian matrix
//
  Epetra_CrsMatrix * HeqProblem::GetMatrix()
  {
    return Matrix;
  }

//
// Returns the Map used to create the solution vector and the Jacobian matrix
//
  Epetra_Map HeqProblem::GetMap()
  {
    return Map;
  }

//
// Calculate and print the 1-norm of the solution vector x along with the parameter value.
// Used to generate the solution graph
//
  void HeqProblem::printSolution(const Epetra_Vector &x, const double conParam)
{
  double n1;                       // temporary variable to hold the value of the norm

  x.Norm1(&n1);                    // calculate the 1-norm of x
  if (outputFilePtr) {             // print out the values of c and ||x||_1
    if (Comm->MyPID()==0) { 
         (*outputFilePtr) << conParam << " " << n1 << endl;
    }
   }
  else 
         cout << "No output file!" << endl;
}

