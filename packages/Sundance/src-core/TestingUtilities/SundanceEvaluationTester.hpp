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

#ifndef SUNDANCE_EVALUATIONTESTER_H
#define SUNDANCE_EVALUATIONTESTER_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "SundanceTestEvalMediator.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceDiscreteFuncElement.hpp"

#include "SundanceEvalManager.hpp"
#include "SundanceOut.hpp"
#include "SundanceExpr.hpp"
#include "SundanceTabs.hpp"
#include "SundanceADCoord.hpp"
#include "SundanceADDerivative.hpp"
#include "SundanceObjectWithVerbosity.hpp"


namespace SundanceTesting
{
  using namespace Sundance;
  using namespace Teuchos;
  using namespace Sundance;
  using namespace Sundance;

  /** 
   *
   */
  class EvaluationTester 
    : public ObjectWithClassVerbosity<EvaluationTester>
  {
  public:
    /** */
    EvaluationTester(const Expr& e, int maxDiffOrder=2);

    /** */
    double evaluate(Array<double>& firstDerivs, 
                    Array<Array<double> >& secondDerivs) const ;
    /** */
    double evaluate(Array<double>& firstDerivs) const ;

    /** */
    double evaluate() const ;
    
    double fdEvaluate(const double& step, const double& tol, 
                      const double& tol2,
                      bool& isOK);

    int numNonzeros() const {return sparsity_->numDerivs();}

    int numNodes() const {return ev_->countNodes();}

  private:
    
    Expr e_;
    RegionQuadCombo rqc_;
    EvalContext context_;
    EvalManager mgr_;
    RCP<AbstractEvalMediator> mediator_;
    mutable TestEvalMediator* tem_;
    const EvaluatableExpr* ev_;
    RCP<SparsitySuperset> sparsity_;
    Map<int, int> unkIDToDiscreteIDMap_;
    int maxDiffOrder_;
  };

}



#endif
