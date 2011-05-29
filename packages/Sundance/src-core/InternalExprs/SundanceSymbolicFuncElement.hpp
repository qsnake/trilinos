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

#ifndef SUNDANCE_SYMBOLICFUNCELEMENT_H
#define SUNDANCE_SYMBOLICFUNCELEMENT_H


#include "SundanceDefs.hpp"
#include "SundanceFuncElementBase.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceSymbolicFuncEvaluator.hpp"
#include "SundanceSymbolicFuncDescriptor.hpp"
#include "SundanceCommonFuncDataStub.hpp"

namespace Sundance
{
  using namespace Sundance;
    class DiscreteFuncElement;
    using namespace Teuchos;

    
    

    /** 
     * SymbolicFuncElement represents a scalar-valued element of a (possibly)
     * list-valued SymbolicFunction. 
     */
    class SymbolicFuncElement : public FuncElementBase,
                                public SymbolicFuncDescriptor,
                                virtual public EvaluatableExpr,
                                public GenericEvaluatorFactory<SymbolicFuncElement, SymbolicFuncElementEvaluator>
    {
    public:
      /** */
      SymbolicFuncElement(const std::string& name, 
        const std::string& suffix,
        const FunctionIdentifier& fid,
        const RCP<const CommonFuncDataStub>& data);
      
      /** virtual destructor */
      virtual ~SymbolicFuncElement() {;}

      /** Append to the set of func IDs present in this expression. */
      void accumulateFuncSet(Set<int>& funcDofIDs, 
        const Set<int>& activeSet) const ;

      /** */
      virtual bool hasTestFunctions() const {return false;}


      /** Specify that expressions involving this function are to be evaluated
       * with this function set to zero. Test functions should always be
       * evaluated at zero. For unknown functions, 
       * substituting zero is appropriate for computing
       * the functional derivatives that arise in a linear problem.
       * */
      void substituteZero() const ;

      /** Specify that expressions involving this function are to be evaluated
       * with this function set to the discrete function (or constant parameter) \f$u_0\f$. 
       * This is appropriate for computing
       * the functional derivatives that arise in a nonlinear expression
       * being linearized about \f$u_0\f$. 
       */
      void substituteFunction(const RCP<DiscreteFuncElement>& u0) const ;

      /** Return the point in function space at which this symbolic 
       * function is to be evaluated. */
      const EvaluatableExpr* evalPt() const {return evalPt_.get();}

      /** Return the point in function space at which this symbolic 
       * function is to be evaluated. */
      EvaluatableExpr* evalPt() {return evalPt_.get();}


      /** */
      bool evalPtIsZero() const ;

      /** */
      const RCP<const CommonFuncDataStub>& commonData() const {return commonData_;}


      /** \name Preprocessing */
      //@{
      /** */
      virtual Set<MultipleDeriv> 
      internalFindW(int order, const EvalContext& context) const ;
          
      /** */
      virtual Set<MultipleDeriv> 
      internalFindV(int order, const EvalContext& context) const ;
          
      /** */
      virtual Set<MultipleDeriv> 
      internalFindC(int order, const EvalContext& context) const ;

      /** */
      virtual RCP<Array<Set<MultipleDeriv> > > 
      internalDetermineR(const EvalContext& context,
                         const Array<Set<MultipleDeriv> >& RInput) const ;
      /** */
      virtual void registerSpatialDerivs(const EvalContext& context, 
                                         const Set<MultiIndex>& miSet) const ;
      //@}
      

      /** Indicate whether the expression is independent of the given 
       * functions */
      virtual bool isIndependentOf(const Expr& u) const ;

      
      /** Indicate whether the expression is linear in the given 
       * functions */
      virtual bool isLinearForm(const Expr& u) const ;
      
      /** */
      virtual RCP<ExprBase> getRcp() {return rcp(this);}
      
    private:
      RCP<const CommonFuncDataStub> commonData_;

      mutable RCP<EvaluatableExpr> evalPt_;

      mutable Array<int> evalPtDerivSetIndices_;
    };
}

#endif
