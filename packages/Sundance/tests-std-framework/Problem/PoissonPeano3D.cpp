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
#include "SundanceEvaluator.hpp"
#include "SundanceLabelCellPredicate.hpp"
#include "SundanceExodusNetCDFMeshReader.hpp"
#include "SundanceElementIntegral.hpp"

/** 
 * Solves the Poisson equation in 3D
 */

#ifdef HAVE_SUNDANCE_PEANO
#ifdef HAVE_SUNDANCE_PEANO_NO_3D

int main(int argc, char** argv)
{
  Sundance::init(&argc, &argv);
  std::cout << "dummy PoissonPeano3D PASSED. Use Peano to run the actual test" << std::endl;
  Sundance::finalize();
  return 0;
}

#else

class PolyFunc : public PointwiseUserDefFunctor0
{
public:
  PolyFunc(int n) : PointwiseUserDefFunctor0("P_" + Teuchos::toString(n), 1, 1), n_(n){}

  /** */
  void eval0(const double* vars, double* f) const ;

private:
  int n_;
};


Expr Poly(int n, const Expr& x)
{
  return  new UserDefOp(x, rcp(new PolyFunc(n)));
}



void PolyFunc::eval0(const double* vars, double* f) const
{
  double y = 1.0;
  double x = vars[0];
  for (int i=0; i<n_; i++)
  {
    double t = 1.0;
    for (int j=0; j<n_; j++) t = t*x;
    y = y + 2.0*t - t - t;
  }
  f[0] = y;
}


CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;})
CELL_PREDICATE(FrontPointTest, {return fabs(x[2]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;})
CELL_PREDICATE(BackPointTest, {return fabs(x[2]-1.0) < 1.0e-10;})

int main(int argc, char** argv)
{
  try
		{
      int depth = 0;
      bool useCCode = false;
      Sundance::ElementIntegral::alwaysUseCofacets() = false;
      Sundance::clp().setOption("depth", &depth, "expression depth");
      Sundance::clp().setOption("C", "symb", &useCCode, "Code type (C or symbolic)");
      Sundance::init(&argc, &argv);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Read the mesh */

      MeshType meshType = new PeanoMeshType3D();
      MeshSource mesher = new PeanoMesher3D(0.0, 0.0,  0.0 , 1.0 , 1.0 , 1.0 , 0.4 , meshType);

      //MeshType meshType = new BasicSimplicialMeshType();
      //MeshSource mesher = new ExodusMeshReader("cube-0.1", meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter faces = new DimensionalCellFilter(2);
      CellFilter side1 = faces.subset(new BottomPointTest());
      CellFilter side2 = faces.subset(new LeftPointTest());
      CellFilter side3 = faces.subset(new FrontPointTest());
      CellFilter side4 = faces.subset(new RightPointTest());
      CellFilter side5 = faces.subset(new BackPointTest());
      CellFilter side6 = faces.subset(new TopPointTest());
      /*
      CellFilter side1 = faces.subset(new LeftPointTest());
      CellFilter side2 = faces.subset(new RightPointTest());
      CellFilter side3 = faces.subset(new TopPointTest());
      CellFilter side4 = faces.subset(new BottomPointTest());
      CellFilter side5 = faces.subset(new FrontPointTest());
      CellFilter side6 = faces.subset(new BackPointTest());
*/
      /* Create unknown and test functions, discretized using second-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(2), "u");
      Expr v = new TestFunction(new Lagrange(2), "v");

      /* Create differential operator and coordinate functions */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr dz = new Derivative(2);
      Expr grad = List(dx, dy, dz);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      Expr z = new CoordExpr(2);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(2);
      QuadratureFamily quad4 = new GaussianQuadrature(4);

      /* Define the weak form */
      //Expr eqn = Integral(interior, (grad*v)*(grad*u) + v, quad);
      
      Expr coeff = 1.0;
#ifdef FOR_TIMING
      if (useCCode)
      {
        coeff = Poly(depth, x);
      }
      else
      {
        for (int i=0; i<depth; i++)
        {
          Expr t = 1.0;
          for (int j=0; j<depth; j++) t = t*x;
          coeff = coeff + 2.0*t - t - t;
        }
      }
#endif
      Expr eqn = Integral(interior, coeff*(grad*v)*(grad*u) +2.0*v, quad2);

      /* Define the Dirichlet BC */
      Expr exactSoln = (x + 1.0)*x - 1.0/4.0;
      Expr h = new CellDiameterExpr();

      Expr bc = EssentialBC( side2 , v*(u-exactSoln), quad4)
         + EssentialBC( side4 , v*(u-exactSoln), quad4 );

      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, v, u, vecType);

#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/aztec-ml.xml"));
#else
      ParameterXMLFileReader reader("aztec-ml.xml");
#endif
      ParameterList solverParams = reader.getParameters();
      std::cerr << "params = " << solverParams << std::endl;


      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      Expr soln = prob.solve(solver);

#ifndef FOR_TIMING

      DiscreteSpace discSpace(mesh, new Lagrange(1), vecType);
      L2Projector proj1(discSpace, exactSoln);
      L2Projector proj2(discSpace, soln-exactSoln);
      L2Projector proj3(discSpace, pow(soln-exactSoln, 2.0));
      Expr exactDisc = proj1.project();
      Expr errorDisc = proj2.project();
//      Expr errorSqDisc = proj3.project();

      std::cerr << "writing fields" << std::endl;
      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("PoissonPeano3d");
      w.addMesh(mesh);
      w.addField("soln", new ExprFieldWrapper(soln[0]));
      w.addField("exact soln", new ExprFieldWrapper(exactDisc));
      w.addField("error", new ExprFieldWrapper(errorDisc));
//      w.addField("errorSq", new ExprFieldWrapper(errorSqDisc));
      w.write();

      std::cerr << "computing error" << std::endl;

      Expr errExpr = Integral(interior, 
                              pow(soln-exactSoln, 2.0),
                              new GaussianQuadrature(4));

      double errorSq = evaluateIntegral(mesh, errExpr);
      std::cerr << "error norm = " << sqrt(errorSq) << std::endl << std::endl;
#else
      double errorSq = 1.0;
#endif
      double tol = 1.0e-10;
      Sundance::passFailTest(sqrt(errorSq), tol);
    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
#endif
#else
int main(int argc, char** argv)
{
  Sundance::init(&argc, &argv);
  std::cout << "dummy PoissonPeano3D PASSED. Use Peano to run the actual test" << std::endl;
  Sundance::finalize();
  return 0;
}

#endif
