#ifndef SUNDANCE_PRODUCTEVALUATOR_H
#define SUNDANCE_PRODUCTEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceBinaryEvaluator.hpp"
#include "SundanceProductExpr.hpp"
#include "Teuchos_TimeMonitor.hpp"


namespace Sundance 
{
/**
 *
 */
class ProductEvaluator : public BinaryEvaluator<ProductExpr>
{
public:
  /** */
  ProductEvaluator(const ProductExpr* expr,
    const EvalContext& context);


  /** */
  virtual ~ProductEvaluator(){;}

  /** */
  virtual void internalEval(const EvalManager& mgr,
    Array<double>& constantResults,
    Array<RCP<EvalVector> >& vectorResults) const ;


  /** */
  TEUCHOS_TIMER(evalTimer, "product evaluation");

private:
  /** */
  enum ProductParity {VecVec, VecConst, ConstVec};
      
  int maxOrder_;
  Array<Array<int> > resultIndex_;
  Array<Array<int> > resultIsConstant_;

  Array<Array<int> > hasWorkspace_;
  Array<Array<int> > workspaceIsLeft_;
  Array<Array<int> > workspaceIndex_;
  Array<Array<int> > workspaceCoeffIndex_;
  Array<Array<int> > workspaceCoeffIsConstant_;

  Array<Array<Array<Array<int> > > > ccTerms_;
  Array<Array<Array<Array<int> > > > cvTerms_;
  Array<Array<Array<Array<int> > > > vcTerms_;
  Array<Array<Array<Array<int> > > > vvTerms_;

  Array<Array<Array<int> > > startingVectors_;
  Array<Array<ProductParity> > startingParities_;
      
      

      
};

}

#endif
