/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
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
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#ifndef SUNDANCE_COLLECTIVEEXCEPTIONCHECK_H
#define SUNDANCE_COLLECTIVEEXCEPTIONCHECK_H

#include "SundanceDefs.hpp"
#include "Teuchos_MPIComm.hpp"


namespace Sundance
{
  using namespace Teuchos;


  /** Call this function upon catching an exception at a point before a collective
   * operation. This function will do an AllReduce in conjunction with calls
   * to either this function or its partner, checkForFailures(), on the
   * other processors. This procedure has the effect of communicating to the other
   * processors that an exception has been detected on this one. */
  void reportFailure(const MPIComm& comm);

  /** Call this function after exception-free completion of a
   * try/catch block preceding a collective operation. This function
   * will do an AllReduce in conjunction with calls to either this
   * function or its partner, reportFailure(), on the other
   * processors. If a failure has been reported by another processor, the
   * call to checkForFailures() will return true and an exception can be thrown. */
  bool checkForFailures(const MPIComm& comm);
}

#endif
