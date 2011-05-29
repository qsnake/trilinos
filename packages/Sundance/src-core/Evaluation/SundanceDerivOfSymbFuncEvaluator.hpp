#ifndef SUNDANCE_DERIVOFSYMBFUNCEVALUATOR_H
#define SUNDANCE_DERIVOFSYMBFUNCEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceUnaryEvaluator.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Sundance 
{
class DerivOfSymbFunc;
class DiscreteFuncElementEvaluator;
    
/**
 *
 */
class DerivOfSymbFuncEvaluator : public UnaryEvaluator<DerivOfSymbFunc>
{
public:
  /** */
  DerivOfSymbFuncEvaluator(const DerivOfSymbFunc* expr,
    const EvalContext& context);

  /** */
  virtual ~DerivOfSymbFuncEvaluator(){;}

  /** */
  virtual void internalEval(const EvalManager& mgr,
    Array<double>& constantResults,
    Array<RCP<EvalVector> >& vectorResults) const ;

  /** We need a specialized resetting method for diff op
   * evaluators that also resets the discrete func evaluators
   * used in the functional chain rule */
  virtual void resetNumCalls() const ;

  /** */
  TEUCHOS_TIMER(evalTimer, "DerivOfSymbFunc evaluation");
private:

  Array<const DiscreteFuncElementEvaluator*> funcEvaluator_;

  int funcMiIndex_;

  bool evalPtIsZero_;

  int constResultIndex_;

  RCP<SparsitySuperset> funcSparsitySuperset_;
}; 
}

#endif
