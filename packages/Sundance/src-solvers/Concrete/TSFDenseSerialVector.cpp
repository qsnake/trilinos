#include "TSFDenseSerialVector.hpp"
#include "Teuchos_Utils.hpp"




using namespace TSFExtended;
using namespace Teuchos;




void DenseSerialVector::setScalar(const double& a)
{
  int n = length();
	double* yy = x();
	for (int i=0; i<n; i++, yy++)
		{
			*yy = a;
		}
}


void DenseSerialVector::negate()
{
  int n = length();
	double* yy = x();
	for (int i=0; i<n; i++, yy++)
		{
			*yy = -(*yy);
		}
}

void DenseSerialVector::add(const DenseSerialVector& other)
{
  int n = length();
  checkLength(other, "add");

	blasObject().AXPY(n, 1.0, other.x(), 1, x(), 1);
}



void DenseSerialVector::subtract(const DenseSerialVector& other)
{
  int n = length();
  checkLength(other, "subtract");

	blasObject().AXPY(n, -1.0, other.x(), 1, x(), 1);
}


void DenseSerialVector::daxpy(const DenseSerialVector& other, const double& a)
{
  int n = length();
  checkLength(other, "daxpy");
	blasObject().AXPY(n, a, other.x(), 1, x(), 1);
}


void DenseSerialVector::eMult(const DenseSerialVector& other)
{
  int n = length();
  checkLength(other, "eMult");
	double* yy = x();
	const double* xx = other.x();	
	for (int i=0; i<n; i++)
		{
			yy[i] *= xx[i];
		}
}


/* added by ptb  */
void DenseSerialVector::abs() 
{
  int n = length();
  double* xx = x();
  for (int i = 0; i < n; i++)
    {
      xx[i] = ::fabs(xx[i]);
    }
}

/* added by ptb  */
double DenseSerialVector::max() const 
{
  int n = length();
  const double* xx = x();
  double ret = xx[0];
  for (int i = 1; i < n; i++)
    {
      if (xx[i] > ret) ret = xx[i];
    }
  return ret;
}

/* added by ptb  */
double DenseSerialVector::min() const 
{
  int n = length();
  const double* xx = x();
  double ret = xx[0];
  for (int i = 1; i < n; i++)
    {
      if (xx[i] < ret) ret = xx[i];
    }
  return ret;
}

/* added by ptb  */
void DenseSerialVector::dotStar(const DenseSerialVector& y, 
                                const DenseSerialVector& z)
{
  int n = length();
  double* xx = x();
  const double* yy = y.x();
  const double* zz = z.x();
  for (int i = 0; i < n; i++)
  {
    xx[i] = yy[i] * zz[i];
  }
}

/* added by ptb  */
void DenseSerialVector::dotSlash(const DenseSerialVector& y, 
                                 const DenseSerialVector& z)
{
  int n = length();
  double* xx = x();
  const double* yy = y.x();
  const double* zz = z.x();
  for (int i = 0; i < n; i++)
  {
    xx[i] = yy[i] / zz[i];
  }
}


void DenseSerialVector::scalarMult(const double& a)
{
  int n = length();
	blasObject().SCAL(n, a, x(), 1);
}

void DenseSerialVector::scalarPow(const double& a)
{
  int n = length();
	if (Teuchos::Utils::chop(a-1.0)==0) return;

	double* yy = x();

	if (Teuchos::Utils::chop(a+1.0)==0)
		{ 
			for (int i=0; i<n; i++)
				{
					yy[i] = 1.0/(yy[i]);
				}
		}
	else if (Teuchos::Utils::chop(a-2.0)==0)
		{ 
			for (int i=0; i<n; i++)
				{
					yy[i] *= yy[i];
				}
		}
	else if (Teuchos::Utils::chop(a-3.0)==0)
		{ 
			for (int i=0; i<n; i++)
				{
					double tmp = yy[i]*yy[i];
					yy[i] *= tmp;
				}
		}
	else if (Teuchos::Utils::chop(a-4.0)==0)
		{ 
			for (int i=0; i<n; i++)
				{
					double tmp = yy[i]*yy[i];
          yy[i] = tmp*tmp;
				}
		}
	else
		{
			for (int i=0; i<n; i++)
				{
					yy[i] = pow(yy[i], a);
				}
		}
}


double DenseSerialVector::dot(const DenseSerialVector& other) const 
{
  int n = length();
  checkLength(other, "dot product");

	double rtn = 0.0;
	const double* xx = x();
	const double* yy = other.x();
	
	rtn = blasObject().DOT(n, xx, 1, yy, 1);
	return rtn;
}

double DenseSerialVector::norm2Squared() const
{
  int n = length();
	double rtn = 0.0;

	rtn = blasObject().NRM2(n, x(), 1);

	return rtn*rtn;
}


double DenseSerialVector::sumElements() const
{
  int n = length();
	double rtn = 0.0;
	const double* xx = x();


	for (int i=0; i<n; i++)
		{
			rtn += xx[i];
		}
	return rtn;
}

double DenseSerialVector::maxNorm() const
{
  int n = length();
	double mx = -1.0e30;
	const double* xx = x();

	for (int i=0; i<n; i++)
		{
      double tmp = fabs(xx[i]);
			if (tmp > mx) mx = tmp;
		}
	return mx;
}



string DenseSerialVector::summary() const
{
  int n = length();
	return "DenseSerialVector[dim=" + Teuchos::toString(n) + "]";
}



			









