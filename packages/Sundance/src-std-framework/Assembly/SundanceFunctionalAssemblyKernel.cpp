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

#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceFunctionalAssemblyKernel.hpp"
#include "SundanceIntegralGroup.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"


using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

using std::setw;
using std::endl;
      
FunctionalAssemblyKernel::FunctionalAssemblyKernel(
  const MPIComm& comm,
  double* value,
  int verb
  )
  : AssemblyKernelBase(verb),
    comm_(comm),
    value_(value),
    localValue_(0.0)
{}

void FunctionalAssemblyKernel::postLoopFinalization()
{
  Tabs tab;
  SUNDANCE_MSG3(verb(), tab << "reducing functional values across processors");
  SUNDANCE_MSG3(verb(), tab << "local value=" << localValue_);
  double globalVal = localValue_;
  
  comm_.allReduce((void*) &localValue_, (void*) &globalVal, 1, 
    MPIComm::DOUBLE, MPIComm::SUM);

  *value_ = globalVal;

  SUNDANCE_MSG3(verb(), tab << "reduced value = " << *value_);
}

void FunctionalAssemblyKernel::fill(
  bool isBC, 
  const IntegralGroup& group,
  const RCP<Array<double> >& localValues) 
{
  Tabs tab;
  SUNDANCE_MSG2(verb(), tab << "adding local increment " << (*localValues)[0]
    << " to local value" << std::endl << tab << "Before: " << localValue_);
  
  localValue_ += (*localValues)[0];
  SUNDANCE_MSG2(verb(), tab << "After: " << localValue_);
}
  

