
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

#include "Teuchos_TimeMonitor.hpp"

using namespace Teuchos;

// Global Timers
RCP<Time> CompTime = TimeMonitor::getNewTimer("Computational Time");
RCP<Time> FactTime = TimeMonitor::getNewTimer("Factorial Time");

// Quadratic function declaration.
double quadFunc( double x );

// Factorial function declaration.
double factFunc( int x );

int main(int argc, char* argv[])
{
  int i;
  double x;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  // Apply the quadratic function.
  for( i=-100; i<100; i++ ) {
    x = quadFunc( (double) i );
  }

  // Apply the factorial function.
  for( i=0; i<100; i++ ) {
    x = factFunc( i );
  }

  // Get a summary from the time monitor.
  TimeMonitor::summarize();

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return 0;
}

/* Evaluate a quadratic function at point x */
double quadFunc( double x )
{
  // Construct a local time monitor, this starts the CompTime timer and will stop when leaving scope.
  Teuchos::TimeMonitor LocalTimer(*CompTime);

  // Evaluate the quadratic function.
  return ( x*x - 1.0 );
}

/* Compute the factorial of x */
double factFunc( int x )
{
  // Construct a local time monitor, this starts the FactTime timer and will stop when leaving scope.
  Teuchos::TimeMonitor LocalTimer(*FactTime);

  // Special returns for specific cases.
  if( x == 0 ) return 0.0;
  if( x == 1 ) return 1.0;

  // Evaluate the factorial function.
  return ( (double) x * factFunc(x-1) );
}

#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure Didasko with:\n"
       "--enable-teuchos");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
#endif
