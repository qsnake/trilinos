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

#ifndef SUNDANCE_EVALUATOR_H
#define SUNDANCE_EVALUATOR_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceSparsitySubset.hpp"
#include "SundanceSet.hpp"
#include "SundanceMultiIndex.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

namespace Sundance 
{
class CoordExpr;

class EvalContext;
  


class EvalManager;

/**
 * Base class for evaluator objects. Each EvaluatableExpr type will 
 * have an associated Evaluator subtype.
 */
class Evaluator : public ObjectWithClassVerbosity<Evaluator>
{
public:
  /** */
  Evaluator();

  /** */
  virtual ~Evaluator(){;}

  /** 
   * Client-level evaluation method. Computes new results on the
   * first call, makes copies on subsequent calls up to the last client, 
   * and finally returns the original result vector upon the 
   * last client's call. 
   */
  void eval(const EvalManager& mgr,
    Array<double>& constantResults,
    Array<RCP<EvalVector> >& vectorResults) const ;

  /** Reset the number of calls to zero. This should be called
   * at the beginning of every new evaluation cycle. */
  virtual void resetNumCalls() const {numCalls_=0;}

  /** */
  virtual void 
  internalEval(const EvalManager& mgr,
    Array<double>& constantResults,
    Array<RCP<EvalVector> >& vectorResults) const = 0 ;

  /** Add one to the number of clients. */
  void addClient() {numClients_++;}

  /** */
  void addConstantIndex(int index, int constantIndex);

  /** */
  void addVectorIndex(int index, int vectorIndex);

      

  /** */
  const Sundance::Map<int, int>& constantIndexMap() const 
    {return constantIndexMap_;}

  /** */
  const Sundance::Map<int, int>& vectorIndexMap() const 
    {return vectorIndexMap_;}
protected:

  /** Return the number of clients that will require results
   * from this evaluator */
  int numClients() const {return numClients_;}

  /** */
  bool isOne(int x) const {return x==1;}

  /** */
  bool isOne(const double& x) const {return isZero(x-1.0);}

  /** */
  bool isZero(const double& x) const {return fabs(x-0.0)<1.0e-15;}

  /** */
  const Array<int>& constantIndices() const {return constantIndices_;}

  /** */
  const Array<int>& vectorIndices() const {return vectorIndices_;}


private:
  int numClients_;

  mutable int numCalls_;

  mutable Array<RCP<EvalVector> > vectorResultCache_;

  mutable Array<double> constantResultCache_;

  Sundance::Map<int, int> constantIndexMap_;

  Sundance::Map<int, int> vectorIndexMap_;

  Array<int> vectorIndices_;

  Array<int> constantIndices_;
};


    

}

#endif
