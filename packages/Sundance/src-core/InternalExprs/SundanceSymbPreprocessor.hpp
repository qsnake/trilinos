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



#ifndef SUNDANCE_SYMBPREPROCESSOR_H
#define SUNDANCE_SYMBPREPROCESSOR_H

#include "SundanceExpr.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceComputationType.hpp"



namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;




class Parameter;


/** */
class SymbPreprocessor 
{
public:
        
  /** */
  static DerivSet setupVariations(const Expr& expr, 
    const Expr& vars,
    const Expr& varEvalPts,
    const Expr& unks,
    const Expr& unkEvalPts,
    const Expr& unkParams,
    const Expr& unkParamEvalPts,
    const Expr& fixedFields,
    const Expr& fixedFieldEvalPts, 
    const Expr& fixedParams,
    const Expr& fixedParamEvalPts, 
    const EvalContext& context,
    const ComputationType& compType);

  /** */
  static DerivSet setupFunctional(const Expr& expr, 
    const Expr& fixedParams,
    const Expr& fixedParamEvalPts,
    const Expr& fixedFields,
    const Expr& fixedFieldEvalPts,
    const EvalContext& context,
  const ComputationType& compType);

  /** */
  static DerivSet setupGradient(const Expr& expr, 
    const Expr& vars,
    const Expr& varEvalPts,
    const Expr& fixedParams,
    const Expr& fixedParamEvalPts,
    const Expr& fixedFields,
    const Expr& fixedFieldEvalPts, 
    const EvalContext& contex,
  const ComputationType& compType);

  /** */
  static DerivSet setupSensitivities(const Expr& expr, 
    const Expr& tests,
    const Expr& unks,
    const Expr& unkEvalPts, 
    const Expr& unkParams,
    const Expr& unkParamEvalPts,
    const Expr& fixedParams,
    const Expr& fixedParamEvalPts,
    const Expr& fixedFields,
    const Expr& fixedFieldEvalPts,
    const EvalContext& context,
  const ComputationType& compType);

  /** */
  static DerivSet setupFwdProblem(const Expr& expr, 
    const Expr& tests,
    const Expr& unks,
    const Expr& unkEvalPts, 
    const Expr& unkParams,
    const Expr& unkParamEvalPts,
    const Expr& fixedParams,
    const Expr& fixedParamEvalPts,
    const Expr& fixedFields,
    const Expr& fixedFieldEvalPts,
    const EvalContext& context,
  const ComputationType& compType);

  /** check the input functions for redundancies and functions
   * of unexpected type. Set evaluation points and collect
   * a set of function IDs. This function is templated so we can use
   * the same code for processing unknown and variational functions */
  template <class T>
  static Set<int> processInputFuncs(const Expr& u, const Expr& u0)
    {
      Set<int> idSet;

      for (int i=0; i<u.size(); i++)
      {
        /* make sure all input functions have the correct type */
        const T* uPtr = dynamic_cast<const T*>(u[i].ptr().get());
        TEST_FOR_EXCEPTION(uPtr==0, RuntimeError, 
          "Unexpected function type error: function " << u[i].toString()
          << " is of type=" << typeid(u).name() 
          << ", but we expected type=" << typeid(T).name());

        /* Add the function's ID to the ID set. While we're here, check
         * to ensure we have no duplicates in the input list. */
        int fid = uPtr->fid().dofID();
        TEST_FOR_EXCEPTION(idSet.contains(fid), RuntimeError,
          "duplicate function in input list " << u.toString());
        idSet.put(fid);


        /* Check that the evaluation point is either a discrete function
         * or a zero expression. */
        RCP<DiscreteFuncElement> u0Ptr
          = rcp_dynamic_cast<DiscreteFuncElement>(u0[i].ptr());
        RCP<ZeroExpr> u0ZeroPtr
          = rcp_dynamic_cast<ZeroExpr>(u0[i].ptr());
        TEST_FOR_EXCEPTION(u0Ptr.get()==NULL && u0ZeroPtr.get()==NULL,
          RuntimeError,
          "evaluation point " << u0[i].toString() << " for func=" << u[i]
          << " is neither a discrete function nor a zero expr");

        /* Set the evaluation point */
        if (u0Ptr.get()==NULL)
        {
          uPtr->substituteZero();
        }
        else
        {
          uPtr->substituteFunction(u0Ptr);
        }
      }

      return idSet;
    }

  /** check the input parameters for redundancies and type. Set 
   * evaluation points and collect
   * a set of parameter IDs. This function is templated so we can use
   * the same code for processing unknown and variational parameters */
  template <class T>
  static Set<int> processInputParams(const Expr& alpha, const Expr& alpha0)
    {
      Set<int> paramID;

      for (int i=0; i<alpha.size(); i++)
      {
        /* ensure everyone has the correct type */
        const T* aPtr = dynamic_cast<const T*>(alpha[i].ptr().get());
        TEST_FOR_EXCEPTION(aPtr==0, RuntimeError,
          "list of purported parameters "
          "contains a function that is not an unknown parameter:"
          << alpha[i].toString());

        int fid = aPtr->fid().dofID();
        TEST_FOR_EXCEPTION(paramID.contains(fid), RuntimeError,
          "duplicate input parameter in list "
          << alpha.toString());
        paramID.put(fid);

        RCP<Parameter> a0Ptr
          = rcp_dynamic_cast<Parameter>(alpha0[i].ptr());
        TEST_FOR_EXCEPTION(a0Ptr.get()==NULL,
          RuntimeError,
          "parameter evaluation point " << alpha0[i].toString()
          << " is not a parameter");
        aPtr->substituteFunction(a0Ptr);
      }
      return paramID;
    }


  /** */
  TEUCHOS_TIMER(preprocTimer, "symbolic preprocessing");
};

/** */
Expr makeZeros(const Expr& e);


}

#endif
