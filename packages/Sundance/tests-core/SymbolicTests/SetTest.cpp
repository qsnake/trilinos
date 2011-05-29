#include "SundanceExpr.hpp"
#include "SundanceStdMathOps.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceSymbolicTransformation.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceParameter.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceDerivSet.hpp"
#include "SundanceRegionQuadCombo.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceStringEvalMediator.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

using Sundance::List;


static Time& totalTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}

int main(int argc, char** argv)
{
  
  try
		{
      GlobalMPISession session(&argc, &argv);

      TimeMonitor t(totalTimer());

      int maxDiffOrder = 0;

      verbosity<SymbolicTransformation>() = 0;
      verbosity<Evaluator>() = 0;
      verbosity<EvalVector>() = 0;
      verbosity<EvaluatableExpr>() = 5;
      Expr::showAllParens() = true;

      EvalVector::shadowOps() = true;

      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);

			Expr u = new UnknownFunctionStub("u");
			Expr w = new UnknownFunctionStub("w");
			Expr v = new TestFunctionStub("v");
			Expr alpha = new UnknownFunctionStub("alpha");

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      Expr u0 = new DiscreteFunctionStub("u0");
      Expr w0 = new DiscreteFunctionStub("w0");
      Expr alpha0 = new DiscreteFunctionStub("alpha0");
      Expr dum;

      Array<Expr> tests;

      tests.append( dx*u );


      for (int i=0; i<tests.length(); i++)
        {
          for (int m=0; m<=maxDiffOrder; m++)
          {
            Out::os() << "=========  diff order "  << m << "========"<<endl;
            RegionQuadCombo rqc(rcp(new CellFilterStub()), 
              rcp(new QuadratureFamilyStub(1)));
            EvalContext context(rqc, m, EvalContext::nextID());
            DerivSet ds = SymbPreprocessor::setupFwdProblem(
              tests[i],
              v, u, u0,
              dum, dum, dum, dum,
              alpha, alpha0,
              context);
            Out::os() << "derivs=" << ds << std::endl;
          }
        }

      TimeMonitor::summarize();
    }
	catch(std::exception& e)
		{
			Out::println(e.what());
		}

  return 0;
  
}
