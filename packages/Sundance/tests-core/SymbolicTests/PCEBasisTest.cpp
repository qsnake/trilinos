#include "SundanceExpr.hpp"
#include "SundanceSpectralBasis.hpp"
#include "SundanceHermiteSpectralBasis.hpp"
#include "Stokhos_HermiteBasis.hpp"
#include "SundanceStokhosBasisWrapper.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"

using std::cout;
using std::exception;
using std::setw;
using namespace Sundance;
using namespace Teuchos;
using namespace Stokhos;



static Time& totalTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}



int main(int argc, char** argv)
{
#ifdef HAVE_SUNDANCE_STOKHOS
  try
  {
    GlobalMPISession session(&argc, &argv);

    TimeMonitor t(totalTimer());

    SpectralBasis h1 = new HermiteSpectralBasis(1, 4);
    SpectralBasis h2 = new Stokhos::HermiteBasis<int, double>(4);

    bool fail = false;
    for (int i=0; i<4; i++)
    {
      for (int j=0; j<4; j++)
      {
        for (int k=0; k<4; k++)
        {
          double c1 = h1.expectation(i,j,k);
          double c2 = h2.expectation(i,j,k);
          double err = fabs(c1 - c2);
          cout << setw(4) << i << setw(4) << j << setw(4) << k
               << setw(16) << c1 << setw(16) << c2 << setw(16) << err;
          if (err > 1.0e-12) 
          {
            cout << " ***** FAILED!" ;
            fail = true;
          }
          cout << std::endl;
        }
      }
    }

    
    TimeMonitor::summarize();
    if (fail) 
    {
      cout << "PCE test FAILED" << std::endl;
      return -1;
    }
    else
    {
      cout << "PCE test PASSED" << std::endl;
    }
  }
	catch(std::exception& e)
  {
    Out::println(e.what());
    return -1;
  }
#else
  std::cout << "test disabled because Stokhos has not been enabled" << std::endl;

#endif
  return 0;  
}

