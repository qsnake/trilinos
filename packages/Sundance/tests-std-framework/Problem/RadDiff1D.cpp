/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "Sundance.hpp"
#include "SundanceProblemTesting.hpp"

int main()
{
  return 0;
}

#ifdef BLARF

/** 
 * Solves the radiation diffusion equation in 1D
 */


double runRadDiff1D(int nx, double* hMean)
{
  LineDomain dom(0.0, 1.0, nx);
  
  VectorType<double> vecType = new EpetraVectorType();

  CellFilter left = dom.left();
  CellFilter right = dom.right();
  CellFilter interior = dom.interior();

  Mesh mesh = dom.mesh();

  int p = 3;
  Expr u = new UnknownFunction(new Lagrange(p), "u");
  Expr v = new TestFunction(new Lagrange(p), "v");

  Expr dx = new Derivative(0);
  Expr x = new CoordExpr(0);

  DiscreteSpace discSpace(mesh, new Lagrange(p), vecType);
  L2Projector projector(discSpace, 1.0+x);
  Expr u0 = projector.project();

  
  QuadratureFamily quad = new GaussianQuadrature(5*p-2);

  Expr eqn = Integral(interior, u*u*u*(dx*v)*(dx*u), quad);

  Expr bc = EssentialBC(left, v*(u-1.0), quad)
    + EssentialBC(right, v*(u-2.0), quad); 

  /* Create a TSF NonlinearOperator object */
  NonlinearProblem prob(mesh, eqn, bc, v, u, u0, vecType);

  
  ParameterXMLFileReader reader("nox-newton.xml");
  ParameterList noxParams = reader.getParameters();
  LinearSolver<double> linSolver = LinearSolverBuilder::createSolver("aztec-ml.xml");
  
  NOXSolver solver(noxParams, linSolver);
  
  // Solve the nonlinear system
  NOX::StatusTest::StatusType status = prob.solve(solver);
  TEST_FOR_EXCEPTION(status != NOX::StatusTest::Converged,
    runtime_error, "solve failed");

  /* check solution */
  Expr exactSoln = pow(15.0*x + 1.0, 0.25);
      
  Expr errExpr = Integral(interior, 
    pow(u0-exactSoln, 2),
    new GaussianQuadrature(4.0*p));
    
  Expr h = new CellDiameterExpr();
  Expr hExpr = Integral(interior, h, new GaussianQuadrature(1));
  Expr AExpr = Integral(interior, 1.0, new GaussianQuadrature(1));
  
  double errorSq = evaluateIntegral(mesh, errExpr);

  double area = evaluateIntegral(mesh, AExpr);
  *hMean = evaluateIntegral(mesh, hExpr)/area;

  return sqrt(errorSq);
}


void fitExp(const Array<double>& h, const Array<double>& err,
  double* p)
{


  Array<double> x(h.size());
  Array<double> y(h.size());
  double xBar = 0.0;
  double yBar = 0.0;
  for (int i=0; i<h.size(); i++)
  {
    x[i] = log(h[i]);
    y[i] = log(err[i]);
    xBar += x[i];
    yBar += y[i];
  }

  xBar /= h.size();
  yBar /= h.size();
  
  double u = 0.0;
  double v = 0.0;
  for (int i=0; i<h.size(); i++)
  {
    u += (x[i]-xBar)*(y[i]-yBar);
    v += pow(x[i]-xBar,2.0);
  }

  *p = u/v;
}


int main(int argc, char** argv)
{
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      Array<int> N;
      Array<double> h;
      Array<double> err;

      for (int i=4; i<=7; i++)
      {
        int nx = lrint(pow(2.0, i));
        double hMean;
        double eps = runRadDiff1D(nx, &hMean);
        N.append(nx*np);
        h.append(hMean);
        err.append(eps);
      }
      for (int i=0; i<N.size(); i++)
      {
        Out::root() << "nx=" << setw(10) << N[i]
                    << " " << setw(20) << setprecision(5) << h[i] 
                    << " " << setw(20) << setprecision(5) << err[i] 
                    << std::endl;
      }
      double p;
      fitExp(h, err, &p);
      Out::root() << "exponent: " << p << std::endl;

    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}


#endif
