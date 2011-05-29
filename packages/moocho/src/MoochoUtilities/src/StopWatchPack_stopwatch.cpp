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

//#include <iostream>
//#include <iosfwd>

#include "StopWatchPack_stopwatch.hpp"
#include "Teuchos_Time.hpp"

double StopWatchPack::seconds(void)
{
  return Teuchos::Time::wallTime();
}

/*

#ifndef _INTEL_CXX

// Implementation using C standard library.
// In MS VC++ 6.0 the precision is only about 0.05 sec.

#include <time.h>

double StopWatchPack::seconds(void)
{
    static const double secs_per_tick = ((double)1.0) / CLOCKS_PER_SEC;
  const clock_t ticks = clock();
  const double sec = ( (double) ticks ) * secs_per_tick;
  //std::cout << "ticks = " << ticks << ", sec = " << sec << std::endl;
    return sec;
}

#else	// _INTEL_CXX implementation

// Windows implementation.

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <assert.h>

namespace {

bool seconds_initialized = false;
LARGE_INTEGER start_count, count_freq;	// counts per sec.

inline void seconds_initialize() {
  if( seconds_initialized ) return;
  // Figure out how often the performance counter increments
  ::QueryPerformanceFrequency( &count_freq );
  // Set this thread's priority as high as reasonably possible to prevent
    // timeslice interruptions
    ::SetThreadPriority( ::GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL );
  // Get the first count.
  TEST_FOR_EXCEPT( !(  QueryPerformanceCounter( &start_count )  ) );
  seconds_initialized = true;
}

}	// end namespace

double StopWatchPack::seconds(void)
{
  seconds_initialize();
  LARGE_INTEGER count;
  QueryPerformanceCounter( &count );
  // "QuadPart" is a 64 bit integer (__int64).  VC++ supports them!
  const double
    sec = (double)( count.QuadPart - start_count.QuadPart ) / count_freq.QuadPart;
  //std::cout << "ticks = " << ticks << ", sec = " << sec << std::endl;
    return sec;
}

#endif	// _INTEL_CXX

*/
