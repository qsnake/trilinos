#include "SundanceExpr.hpp"
#include "SundanceStdMathOps.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceSymbolicTransformation.hpp"
#include "SundancePointwiseUserDefFunctor.hpp"
#include "SundanceUserDefOp.hpp"
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
#include "SundanceEvaluationTester.hpp"

using namespace Sundance;
using namespace SundanceTesting;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;


static Time& totalTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}


class MyFunc : public PointwiseUserDefFunctor2
{
public:
  MyFunc(int n) : PointwiseUserDefFunctor2("MyFunc", n, 1){;}
  virtual ~MyFunc(){;}
  void eval1(const double* vars, double* f, double* df) const ;
  void eval0(const double* vars, double* f) const ;
  virtual void eval2(const double* in, double* outVals, double* outDerivs,
    double* outDerivs2) const ;
private: 
  int n_;
};



void MyFunc::eval0(const double* vars, double* f) const
{
  f[0] = 0.0;
  for (int i=0; i<n_; i++)
  {
    for (int j=0; j<n_; j++)
    {
      f[0] += (i+1)*(j+1)*vars[i] * vars[j];
    }
  }
}


void MyFunc::eval1(const double* vars, double* f, double* df) const
{
  f[0] = 0.0;
  for (int i=0; i<n_; i++)
  {
    df[i] = 0.0;
    for (int j=0; j<n_; j++)
    {
      df[i] += (i+1)*(j+1)*vars[j];
      f[0] += (i+1)*(j+1)*vars[i] * vars[j];
    }
  }
}
void MyFunc::eval2(const double* vars, double* f, double* df, double* outDerivs2) const
{
  f[0] = 0.0;
  for (int i=0; i<n_; i++)
  {
    df[i] = 0.0;
    for (int j=0; j<n_; j++)
    {
      df[i] += (i+1)*(j+1)*vars[j];
      f[0] += (i+1)*(j+1)*vars[i] * vars[j];
      outDerivs2[i*n_+j] = (i+1)*(j+1);
    }
  }
}




int main(int argc, char** argv)
{
  
  try
		{
      GlobalMPISession session(&argc, &argv);
      Tabs tabs;
      TimeMonitor timer(totalTimer());

      //       verbosity<SymbolicTransformation>() = 0;
      //       verbosity<EvaluationTester>() = 5;
      //       verbosity<Evaluator>() = 5;
      //       verbosity<EvalVector>() = 5;
      //       verbosity<EvaluatableExpr>() = 5;
      //       verbosity<AbstractEvalMediator>() = 5;
      Expr::showAllParens() = true;

      EvalVector::shadowOps() = false;

      Time stopwatch("test");

      for (int n=1; n<0; n++)
        {
          double t0 = stopwatch.wallTime();
          stopwatch.start();
          int depth = n;

          ADField U(ADBasis(1), sqrt(2.0));
          Expr u = new TestUnknownFunction(U, "u");
          Expr x = new CoordExpr(0);

          Expr expr = u + x;

          for (int i=0; i<depth; i++)
            {
              expr = expr * expr;
            }

          EvaluationTester tester(expr);
          Array<double> df;
          Array<Array<double> > df2;
          tester.evaluate(df, df2);
          stopwatch.stop();
          double t1 = stopwatch.wallTime();
          //         std::cerr << n << "   "  << tester.numNodes() << "    " << t1-t0 << std::endl;

        }

      for (int n=2; n<54; n+=4)
        {


          int depth = n;

          ADField U(ADBasis(1), sqrt(2.0));
          Array<Expr> unks(n);
          Array<ADField> Unks(n);
          for (int i=0; i<n; i++) 
            {
              Unks[i] = ADField(ADBasis(1), sqrt(double(1+i)));
              unks[i] = new TestUnknownFunction(Unks[i], "u" + Teuchos::toString(i));
            }
//          Expr uList = new ListExpr(unks);
//          Expr userDefExpr = new UserDefOp(uList, rcp(new MyFunc(n)));

          Expr x = new CoordExpr(0);

          stopwatch.start();

          Expr expr = x;
          for (int i=0; i<depth; i++)
            {
              for (int j=0; j<depth; j++)
                {
                  if (rand() < RAND_MAX/2)
                  expr = expr + unks[i]*unks[j];
                }
            }

          double t0 = stopwatch.wallTime();
          EvaluationTester tester(expr);
          stopwatch.stop();
          double t1 = stopwatch.wallTime();
          Array<double> df;
          Array<Array<double> > df2;


          tester.evaluate(df, df2);


//          stopwatch.start();
//          EvaluationTester udTester(userDefExpr);
//          stopwatch.stop();
//          double t2 = stopwatch.wallTime();
          std::cerr << n << "   "  
               << tester.numNodes() << "    " << t1-t0 << std::endl;

        }
      TimeMonitor::summarize();
    }
	catch(std::exception& e)
		{
			Out::println(e.what());
		}


  
}
