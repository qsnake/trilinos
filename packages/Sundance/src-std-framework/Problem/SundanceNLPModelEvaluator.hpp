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
#ifndef SUNDANCE_SUNDANCENLPMODELEVALUATOR_H
#define SUNDANCE_SUNDANCENLPMODELEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceDefs.hpp"

#ifdef HAVE_SUNDANCE_MOOCHO

#include "SundanceNonlinearProblem.hpp"
#include "SundanceLinearProblem.hpp"
#include "SundanceFunctional.hpp"
#include "SundanceModelEvaluatorBase.hpp"


namespace Thyra
{
  class SundanceNLPModelEvaluator : public SundanceModelEvaluator
  {
  public:
    /** */
    SundanceNLPModelEvaluator(const VectorType<double>& vecType);

    /** */
    void initialize(const Expr& paramExpr,
                    const Expr& stateExpr,
                    const Expr& stateExprVal,
                    const NonlinearProblem& prob, 
                    const Array<LinearProblem>& sensProb,
                    const Functional& objective);


    /** */
    virtual Vector<double> getInitialState() const 
    {
      Vector<double> rtn = stateSpace().createMember();
      rtn.setToConstant(0.0);
      return rtn;
    }

    
    /** */
    virtual Vector<double> getInitialParameters() const ;


    /** */
    void setInitialParameters(const Array<double>& a);

    /** */
    void internalEvalModel(const Vector<double>& stateVec,
                           const Vector<double>& params,
                           Vector<double>& resid,
                           double& objFuncVal,
                           LinearOperator<double>& df_dx,
                           Array<Vector<double> >& df_dp,
                           Vector<double>& dg_dp_T,
                           Vector<double>& dg_dx_T) const ;
    

    /** */
    VectorSpace<double> paramSpace() const 
    {
      return paramSpace_;
    }        
           
    /** */
    VectorSpace<double> stateSpace() const 
    {
      VectorSpace<double> rtn = createW().domain();
      return rtn;
    }
           
    /** */
    VectorSpace<double> constraintSpace() const 
    {
      VectorSpace<double> rtn = createW().range();
      return rtn;
    }
           
    /** */
    LinearOperator<double> createW() const 
    {
      static LinearOperator<double> J = prob_.allocateJacobian();
      TEST_FOR_EXCEPTION(J.ptr().get()==0, RuntimeError,
                         "null Jacobian");
      return J;
    }

    /** */
    Array<double> parameters() const ;

    /** */
    Expr stateVariable() const {return stateExpr_;}


    /** */
    Expr solveForward(const ParameterList& fwdParams) const ;


    /** */
    Array<double> paramArray(const ParameterList& params,
                             const std::string& paramName) const ;

    /** */
    void setContinuationParameters(const Expr& contParams) {contParams_ = contParams;}

    /** */
    void setFinalContinuationValues(const Expr& finalContParams) 
    {finalContParams_ = finalContParams;}

    /** */
    Expr continuationParameters(int i) const {return contParams_[i];}
    /** */
    int numContinuationParameters() const {return contParams_.size();}

    /** */
    Expr finalContinuationValues(int i) const {return finalContParams_[i];}
    

  private:
    VectorSpace<double> paramSpace_;

    Array<double> initParams_;

    mutable Expr paramExpr_;

    mutable Expr stateExpr_;

    mutable NonlinearProblem prob_;

    Array<LinearProblem> sensProb_;
    
    Functional obj_;
    
    FunctionalEvaluator objEval_;

    Expr contParams_;

    Expr finalContParams_;


  };
}


#endif
#endif
