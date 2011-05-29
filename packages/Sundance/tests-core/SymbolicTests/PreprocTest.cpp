#include "SundanceSymbPreprocessor.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "Teuchos_TestingHelpers.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using std::cout;
using std::exception;

using Sundance::List;


#define TEST_THROW(code, passFail) \
  TEUCHOS_TEST_THROW( code, std::exception, Out::os(), passFail)

#define TEST_NOTHROW(code, passFail) \
  TEUCHOS_TEST_NOTHROW( code, Out::os(), passFail)

bool validateFuncTypeChecking()
{
  Expr ux = new UnknownFunctionStub("ux");
  Expr vx = new TestFunctionStub("vx");
  Expr uy = new UnknownFunctionStub("uy");
  Expr vy = new TestFunctionStub("vy");
  Expr uz = new UnknownFunctionStub("uz");
  Expr vz = new TestFunctionStub("vz");

  Expr v = List(vx,vy,vz);
  Expr u = List(ux,uy,uz);

  Expr mixup = List(vx, uy, vz); // mix of test & unknown
  Expr dup = List(vx, vx, vz); // list with duplicates

  bool passFail = true;
  /* */
  Out::os() << "Testing detection of mixed-up function types" << std::endl;
  TEST_THROW(
    SymbPreprocessor::processInputFuncs<UnknownFuncElement>(mixup, makeZeros(v)), passFail);
  
  /* */
  Out::os() << "Testing detection of duplicated functions" << std::endl;
  TEST_THROW(
    SymbPreprocessor::processInputFuncs<SymbolicFuncElement>(dup, makeZeros(v)), passFail);
  
  /* */
  Out::os() << "Testing detection of invalid evaluation points" << std::endl;
  TEST_THROW(
    SymbPreprocessor::processInputFuncs<SymbolicFuncElement>(v, u), passFail);
  

  /* */
  Out::os() << "Testing processing of good input" << std::endl;
  TEST_NOTHROW(
    SymbPreprocessor::processInputFuncs<SymbolicFuncElement>(v, makeZeros(v)), passFail);
  
  return passFail;
} 

int main(int argc, char** argv)
{
  bool pass = true;
  
  try
		{
      GlobalMPISession session(&argc, &argv);

      pass = pass && validateFuncTypeChecking();
    }
	catch(std::exception& e)
		{
      pass = false;
			Out::println(e.what());
		}

  if (pass)
  {
    Out::os() << "test PASSED" << std::endl;
    return 0;
  }
  else 
  {
    Out::os() << "test FAILED" << std::endl;
    return -1;
  }

  
}
