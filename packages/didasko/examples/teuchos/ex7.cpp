
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

#include "Didasko_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#endif
#if defined(HAVE_DIDASKO_TEUCHOS)

#include "Teuchos_CommandLineProcessor.hpp"

int main(int argc, char* argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Creating an empty command line processor looks like:
  Teuchos::CommandLineProcessor My_CLP;

  /* To set and option, it must be given a name and default value.  Additionally,
     each option can be given a help string.  Although it is not necessary, a help
     string aids a users comprehension of the acceptable command line arguments.
     Some examples of setting command line options are:
  */
  // Set an integer command line option.
  int NumIters = 1550;
  My_CLP.setOption("iterations", &NumIters, "Number of iterations");
  // Set a double-precision command line option.
  double Tolerance = 1e-10;    
  My_CLP.setOption("tolerance", &Tolerance, "Tolerance");
  // Set a string command line option.
  string Solver = "GMRES";
  My_CLP.setOption("solver", &Solver, "Linear solver");
  // Set a boolean command line option.    
  bool Precondition;
  My_CLP.setOption("precondition","no-precondition",
		   &Precondition,"Preconditioning flag");

  /* There are also two methods that control the strictness of the command line processor.
     For a command line processor to be sensitive to any bad command line option that it 
     does not recognize use:
  */
  My_CLP.recogniseAllOptions(false);
  
  /* Then, if the parser finds a command line option it doesn't recognize, it will
     throw an exception.  To prevent a command line processor from throwing an exception 
     when it encounters a unrecognized option or help is printed, use:
  */
  My_CLP.throwExceptions(false);
  
  //Finally, to parse the command line, argc and argv are passed to the parse method:
  My_CLP.parse( argc, argv );

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure Didasko with:");
  puts("--enable-teuchos");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
#endif
