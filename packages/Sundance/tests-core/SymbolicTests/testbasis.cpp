#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "Sundance.hpp"
#include "SundanceSpectralBasis.hpp"
#include "SundanceSpectralExpr.hpp"

int main(int argc, void** argv)
{
  try
    {
      Sundance::init(&argc, &argv);
      
      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      int ndim = 2;
      int order = 2;

      SpectralBasis SB(ndim, order);
     
      Array<Expr> Ex(6);
      Ex[0] = x;
      Ex[1] = y;
      Ex[2] = x*y;
      Ex[3] = 0.0;
      Ex[4] = x*x;
      Ex[5] = y*y;

      Expr SE = new SpectralExpr(SB, Ex);

      const SpectralExpr* se = dynamic_cast<const SpectralExpr*>(SE.ptr().get());
      SpectralBasis basis = se->getSpectralBasis();

      double e;

      for (int i=0; i< basis.nterms(); i++)
	for (int j=0; j< basis.nterms(); j++)
	  for(int k=0; k< basis.nterms(); k++)
	    {
	      e = basis.expectation(basis.getElement(i), basis.getElement(j), basis.getElement(k));
	      cout << i << " " << j  << " " << k << " " << e << std::endl;
	    }

    }
  
  catch(std::exception& e)
    {
      Sundance::handleException(e);
    }
  Sundance::finalize();
  return Sundance::testStatus();
}

