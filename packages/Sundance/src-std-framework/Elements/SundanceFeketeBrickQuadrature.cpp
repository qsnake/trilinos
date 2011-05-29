#include "SundanceFeketeBrickQuadrature.hpp"
#include "SundanceOut.hpp"
#include "SundanceGaussLobatto1D.hpp"
#include "SundanceTabs.hpp"

using namespace Sundance;
using namespace Teuchos;

void FeketeBrickQuadrature::getPoints(int order, Array<double>& wgt, Array<
		double>& x, Array<double>& y, Array<double>& z)
{
	int p = order + 3;
	p = p + (p % 2);
	int nNodes = p / 2;
	GaussLobatto1D rule(nNodes, 0.0, 1.0);
	Array<double> d1 = rule.nodes();
	Array<double> d2 = d1;
	Array<double> d3 = d1;
	Array<double> w = rule.weights();
	int n = rule.nPoints();

	wgt.resize(n * n * n);
	x.resize(n * n * n);
	y.resize(n * n * n);
	z.resize(n * n * n);

	int k = 0;
	for (int i = 0; i < n; i++)
	{
		double p = d1[i];
		for (int j = 0; j < n; j++)
		{
			double q = d2[j]; //similar to the p value, caz we have quad
			for (int l = 0; l < n; l++, k++)
			{
				double r = d3[l];
				x[k] = p;
				y[k] = q;
				z[k] = r;
				wgt[k] = w[i] * w[j] * w[l];
			}
		}
	}
}

bool FeketeBrickQuadrature::test(int p)
{
	Array<double> w;
	Array<double> x;
	Array<double> y;
	Array<double> z;

	getPoints(p, w, x, y, z);
	bool pass = true;
	for (int a = 0; a <= p; a++)
	{
		int bMax = p - a;
		for (int b = 0; b <= bMax; b++)
		{
			int cMax = bMax - b;
			for (int c = 0; c <= cMax; c++)
			{
				double sum = 0.0;
				for (int q = 0; q < w.length(); q++)
				{
					sum += w[q] * pow(x[q], (double) a) * pow(y[q], (double) b)
							* pow(z[q], (double) c);
				}
				double err = fabs(sum - exact(a, b, c));
				bool localPass = err < 1.0e-14;
				pass = pass && localPass;
				if (!localPass)
				{
					fprintf(
							stderr,
							"order=%d m (%d, %d, %d) q=%22.15g exact=%22.15g\n",
							p, a, b, c, sum, exact(a, b, c));
					std::cerr << "error = " << err << std::endl;
				}
			}
		}
	}
	return pass;
}

double FeketeBrickQuadrature::exact(int a, int b, int c)
{
	return 1.0 / (a + 1) / (b + 1) / (c + 1);
}

