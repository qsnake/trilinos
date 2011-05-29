#include "SundanceExpr.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnknownFunctionStub.hpp"
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
#include "SundanceFunctionalPolynomial.hpp"

#include "SundanceSymbolicTransformation.hpp"
#include "SundanceProductTransformation.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceEquationSet.hpp"

using Sundance::List;
using namespace Sundance;
using namespace Teuchos;


RCP<ScalarExpr> expr2scalar(const Expr& e)
{
  return rcp_dynamic_cast<ScalarExpr>(e[0].ptr());
}

RCP<FunctionalPolynomial> expr2poly(const Expr& e)
{
  return FunctionalPolynomial::toPoly(expr2scalar(e));
}

int main(int argc, char** argv)
{
  
  try
		{
      GlobalMPISession session(&argc, &argv);

      Expr::showAllParens() = true; 
      ProductTransformation::optimizeFunctionDiffOps() = true;

      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);

      Expr u = new UnknownFunctionStub("u");
      Expr v = new UnknownFunctionStub("v");
      Expr w = new UnknownFunctionStub("w");
      
      RCP<FunctionalPolynomial> p 
        = rcp(new FunctionalPolynomial(expr2scalar(u)));
      
      RCP<FunctionalPolynomial> q 
        = rcp(new FunctionalPolynomial(expr2scalar(w)));
      
      RCP<FunctionalPolynomial> r
        = rcp(new FunctionalPolynomial(expr2scalar(w)));

      p = p->multiplyPoly(expr2poly(v).get());
      p = p->multiplyPoly(expr2poly(v).get());
      p = p->multiplyPoly(expr2poly(u).get());

      q = q->multiplyPoly(expr2poly(u).get());
      q = q->multiplyPoly(expr2poly(v).get());

      p = p->addPoly(q.get(), 1);

      p = p->addPoly(r.get(), 1);

          
      Out::os() << "p = " << p->toXML() << std::endl;
      p->toText(Out::os(), false) << std::endl;

      TimeMonitor::summarize();      
    }
	catch(std::exception& e)
		{
      Out::os() << e.what() << std::endl;
		}

  
}
