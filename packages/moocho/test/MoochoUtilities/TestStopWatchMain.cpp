// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "StopWatchPack_stopwatch.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

int main(int argc, char* argv[])
{

  using std::cout;
  using std::endl;
  using StopWatchPack::stopwatch;
  
  bool success = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  try {
    
    stopwatch timer;
    
    // Finding min resolution.
    double min_resolution = 0.0; // in seconds
    int total_num_calls   = 0;
    {
      cout	<< "\n*** Measuring miminum resolution.\n";
      timer.start();
      double last_time = timer.read();
      const int max_num_samples = 20;
      int num_samples = 0;
      int num_calls = 0;
      while( num_samples < max_num_samples ) {
        double time = timer.read();
        num_calls++;
        if( time - last_time > 0.0 ) {
          cout	<< "time_diff = " << time - last_time
                << ", num_calls = " << num_calls << endl;
          min_resolution += time - last_time;
          ++total_num_calls;
          last_time = time;
          num_calls = 0;
          num_samples++;
        }
      }
      min_resolution /= total_num_calls;
    }
    
    std::cerr << "Minimum stopwatch resolution = " << min_resolution << " sec\n";
    
    // Finding increasing resolution.
    {
      cout	<< "\n*** Measuring increasing resolution.\n";
      timer.start();
      double start_time = timer.read(), last_time = start_time;
      const int max_num_samples = 20;
      int num_samples = 0;
      int num_calls = 0;
      while( num_samples < max_num_samples ) {
        double time = timer.read();
        num_calls++;
        if( time - last_time > 0.0 ) {
          cout	<< "time = " << time - start_time
                << ", num_calls = " << num_calls << endl;
          last_time = time;
          num_calls = 0;
          num_samples++;
        }
      }
    }

  } // end try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  
  if(success)
    cout << "\nEnd Result: TEST PASSED" << std::endl;
  
  return 0;
  
}
