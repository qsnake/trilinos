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

#ifndef STPWATCH_H
#define STPWATCH_H

#include "Moocho_ConfigDefs.hpp"

namespace StopWatchPack {

/** @name namespace StopWatchPack
  *
  * @memo Package for CPU timing.
  */
//@{

double seconds(void);

/** \brief Simple stopwatch object.
  */
class stopwatch {
public:

  /// Initializes of not running.
  stopwatch() : running_(false), last_time_(0.0), total_(0.0)
  {}

  /// Returns true if <tt>this</tt> is currently timming.
  bool is_running() const {
    return running_;
  }

  /// Starts timing if it has already not been started.
  void start() {
    if (!running_) {
      last_time_ = seconds();
      running_ = true;
    }
  }

  /// Stops timing and returns the time (sec.) since start() was called
  double stop()  {
    if (running_) {
      total_ += seconds() - last_time_; 
      running_ = false;
    }
    //std::cout << "total = " << total_ << std::endl;
    return total_; 
  }

  /// Stops and resets the clock if it is running.
  void reset() {
    running_ = false;
    last_time_ = 0.0;
    total_ = 0.0;
  }

  /// Reads the elapsed time (sec.) and leaves the clock running.
  double read() {
    if (running_) {
      double curr_time = seconds();
      total_ += curr_time - last_time_;
      last_time_ = curr_time;
    }
    return total_;
  }
  
private:
  bool running_;
  double last_time_;
  double total_;
};

//	end namespace StopWatchPack 
//@}

}  // end namespace StopWatchPack 

#endif // STPWATCH_H
