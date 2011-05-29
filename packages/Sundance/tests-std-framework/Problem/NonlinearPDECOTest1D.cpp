#include "Sundance.hpp"

#include "SundanceZeroExpr.hpp"

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;}) 
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;}) 


int main(int argc, char** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();
      int npx = -1;
      int npy = -1;
      PartitionedRectangleMesher::balanceXY(np, &npx, &npy);
      TEST_FOR_EXCEPT(npx < 1);
      TEST_FOR_EXCEPT(npy < 1);
      TEST_FOR_EXCEPT(npx * npy != np);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      MeshType meshType = new BasicSimplicialMeshType();
      int nx = 128;
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx, meshType); 
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter(); 
      CellFilter pts = new DimensionalCellFilter(0);
      
      CellFilter left = pts.subset(new LeftPointTest());
      CellFilter right = pts.subset(new RightPointTest());

      /* Create unknown function */
      Expr u = new UnknownFunction(new Lagrange(1), "u");
      Expr lambda = new UnknownFunction(new Lagrange(1), "lambda");
      Expr alpha = new UnknownFunction(new Lagrange(1), "alpha");
      
      Expr dx = new Derivative(0);
      Expr grad = dx;
      
      Expr x = new CoordExpr(0);
      
      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(4);

      WatchFlag watch("watch me");
      watch.setParam("integration setup", 6);
      watch.deactivate();

      const double pi = 4.0*atan(1.0);

      double R = 0.001;

      /* uStar, spliced in from Mathematica */
      Expr uStar = 2*pow(pi,2)*R*
        pow(cos(pi*x),2) + 
        sin(pi*x) - 
        pow(pi,2)*R*
        pow(sin(pi*x),2) + 
        2*R*pow(sin(pi*x),3);
      
      Expr mismatch = u-uStar;
      Expr fit = Integral(interior, 0.5*mismatch*mismatch, quad, watch);
      Expr reg = Integral(interior, 0.5*R*alpha*alpha, quad, watch);

      Expr g = pi*pi*u + u*u;

      Expr constraint = Integral(interior, (grad*u)*(grad*lambda) - lambda*g + lambda*alpha, quad, watch);

      Expr lagrangian = fit + reg + constraint;

      Expr bc = EssentialBC(left+right, 
        lambda*u, quad, watch);

      BasisFamily L1=new Lagrange(1);
      DiscreteSpace discSpace(mesh, List(L1, L1, L1), vecType);
      Expr W0 = new DiscreteFunction(discSpace, 1.0);
      Expr u0 = W0[0];
      Expr lambda0 = W0[1];
      Expr alpha0 = W0[2];
      Expr W = List(u, lambda, alpha);

      Expr dummy;

      Functional L(mesh, lagrangian, bc, vecType);
      
      NonlinearProblem prob 
        = L.nonlinearVariationalProb(W, W0, W, W0, dummy, dummy);

      ParameterXMLFileReader reader("nox-amesos.xml");

      ParameterList noxParams = reader.getParameters();

      NOXSolver solver(noxParams);
      
      NOX::StatusTest::StatusType status = prob.solve(solver);
      TEST_FOR_EXCEPTION(status != NOX::StatusTest::Converged,
        runtime_error, "solve failed");


      /* Write the field in VTK format */
      FieldWriter w = new MatlabWriter("pdeco1D");
      w.addMesh(mesh);
      w.addField("u", new ExprFieldWrapper(W0[0]));
      w.addField("lambda", new ExprFieldWrapper(W0[1]));
      w.addField("alpha", new ExprFieldWrapper(W0[2]));
      w.write();



      Expr uExact = sin(pi*x);
      Expr alphaExact = uExact*uExact;
      Expr lambdaExact = -R*alphaExact;

      Expr uErrExpr = Integral(interior, 
        pow(u0-uExact, 2),
        new GaussianQuadrature(8));

      Expr alphaErrExpr = Integral(interior, 
        pow(alpha0-alphaExact, 2),
        new GaussianQuadrature(8));

      Expr lambdaErrExpr = Integral(interior, 
        pow(lambda0-lambdaExact, 2),
        new GaussianQuadrature(8));

      
      double uErrorSq = evaluateIntegral(mesh, uErrExpr);
      std::cerr << "u error norm = " << sqrt(uErrorSq) << std::endl << std::endl;

      double alphaErrorSq = evaluateIntegral(mesh, alphaErrExpr);
      std::cerr << "alpha error norm = " << sqrt(alphaErrorSq) << std::endl << std::endl;

      double lambdaErrorSq = evaluateIntegral(mesh, lambdaErrExpr);
      std::cerr << "lambda error norm = " << sqrt(lambdaErrorSq) << std::endl << std::endl;

      double err = sqrt(uErrorSq + lambdaErrorSq + alphaErrorSq);

      double tol = 0.01;
      Sundance::passFailTest(err, tol);
    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}
  Sundance::finalize();
  return Sundance::testStatus(); 
}
