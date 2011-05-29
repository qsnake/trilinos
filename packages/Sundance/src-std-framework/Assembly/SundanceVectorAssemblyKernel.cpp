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
#include "SundanceVectorAssemblyKernel.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"


using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using namespace TSFExtended;
using std::setw;
using std::endl;
      
VectorAssemblyKernel::VectorAssemblyKernel(
  const Array<RCP<DOFMapBase> >& dofMap,
  const Array<RCP<Array<int> > >& isBCIndex,
  const Array<int>& lowestLocalIndex,
  Array<Vector<double> >& b,
  bool partitionBCs,
  int verb
  )
  : VectorFillingAssemblyKernel(dofMap, isBCIndex, lowestLocalIndex,
    b, partitionBCs, verb)
{}


void VectorAssemblyKernel::fill(
  bool isBC, 
  const IntegralGroup& group,
  const RCP<Array<double> >& localValues) 
{
  Tabs tab0;
  SUNDANCE_MSG1(verb(), tab0 << "in VectorAssemblyKernel::fill()");

  TEST_FOR_EXCEPT(!group.isOneForm());

  bool useCofacets = group.usesMaximalCofacets();

  if (group.isOneForm())
  {
    insertLocalVectorBatch(isBC, useCofacets, 
      group.testID(), group.testBlock(), group.mvIndices(),
      *localValues);
  }

  SUNDANCE_MSG1(verb(), tab0 << "done VectorAssemblyKernel::fill()");
}
  

void VectorAssemblyKernel:: prepareForWorkSet(
  const Array<Set<int> >& requiredTests,
  const Array<Set<int> >& /* requiredUnks */,
  RCP<StdFwkEvalMediator> mediator)
{
  Tabs tab0;
  SUNDANCE_MSG1(verb(), tab0 
    << "in VectorAssemblyKernel::prepareForWorkSet()");
  IntegrationCellSpecifier intCellSpec = mediator->integrationCellSpec();
  buildLocalDOFMaps(mediator, intCellSpec, requiredTests);
  SUNDANCE_MSG1(verb(), tab0 << "done VectorAssemblyKernel::prepareForWorkSet()");
}
