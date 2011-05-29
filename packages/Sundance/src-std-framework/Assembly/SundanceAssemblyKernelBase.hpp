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

#ifndef SUNDANCE_ASSEMBLYKERNELBASE_H
#define SUNDANCE_ASSEMBLYKERNELBASE_H

#include "SundanceDefs.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceIntegralGroup.hpp"
#include "SundanceElementIntegral.hpp"
#include "SundanceGrouperBase.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceStdFwkEvalMediator.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "TSFLoadableVector.hpp"
#include "TSFLoadableMatrix.hpp"
#include "TSFLinearOperatorDecl.hpp"
#include "TSFVectorDecl.hpp"
#include "TSFVectorType.hpp"
#include "Teuchos_HashSet.hpp"
#include "Teuchos_ParameterList.hpp"
#include "TSFIncrementallyConfigurableMatrixFactory.hpp"
#include "TSFCollectivelyConfigurableMatrixFactory.hpp"
#include "TSFPartitionedMatrixFactory.hpp"
#include "TSFPartitionedToMonolithicConverter.hpp"

namespace Sundance
{
using namespace Teuchos;

class StdFwkEvalMediator;
class IntegralGroup;

/** 
 * AssemblyKernelBase abstracts the operations that must be done in an
 * assembly loop. Regardless of whether the assembly loop is doing
 * matrix/vector fill, vector fill, functional/gradient evaluation,
 * or functional evaluation, the assembly loop will involve
 * <ul>
 * <li> A preprocessing step before the main assembly loop
 * <li> A preprocessing step before each work set is started
 * <li> A fill step after each integral group is done
 * <li> A postprocessing step after the main assembly loop is done
 * </ul>
 * The first of these is done by the subclass constructor. The others
 * are done using the pure virtual functions of this class.
 *
 * It is assumed that any data structures to be filled -- such as a matrix,
 * a vector, or simply a number -- are stored internally in the assembly
 * kernel subclass, and that they persist between preprocessing and fill 
 * calls.
 */
class AssemblyKernelBase
{
public:
  /** */
  AssemblyKernelBase(int verb) : verb_(verb) {;}

  /** */
  virtual ~AssemblyKernelBase(){;}

  /**  
   * Do preprocessing steps needed before integrating the current
   * work set. 
   *
   * The default implementation does nothing. 
   */
  virtual void prepareForWorkSet(
    const Array<Set<int> >& requiredTests,
    const Array<Set<int> >& requiredUnks,
    RCP<StdFwkEvalMediator> mediator) {;}

  /** 
   * Adds the results of the current integral
   * group into the assembly results. 
   * \param isBC whether the current group is a replace-style 
   * boundary condition
   * \param group the current integral group
   * \param localValues the results of integrating the current integral group
   */
  virtual void fill(bool isBC,
    const IntegralGroup& group,
    const RCP<Array<double> >& localValues) = 0 ;  

  /** 
   * Hook to do any finalization steps after the main assembly loop, 
   * for example, doing an all-reduce on locally computed functional values. 
   * The default implementation does nothing. */
  virtual void postLoopFinalization() {;}

  /** verbosity level */
  int verb() const {return verb_;}

  /** set verbosity level.
   * (This function needs to be virtual because certain subclasses need specialized
   * implementations that propagate verbosity to children 
  */
  virtual void setVerbosity(int verb) {verb_=verb;}

private:
  int verb_;
};



}




#endif
