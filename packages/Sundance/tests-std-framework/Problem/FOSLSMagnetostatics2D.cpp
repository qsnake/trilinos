#include "Sundance.hpp"

#include "SundanceZeroExpr.hpp"
CELL_PREDICATE(TopLeftTest, {return fabs(x[0]) < 1.0e-10 && x[1]>=0.5;}) 
CELL_PREDICATE(BottomLeftTest, {return fabs(x[0]) < 1.0e-10 && x[1]<=0.5;}) 
CELL_PREDICATE(TopRightTest, {return fabs(x[0]-1.0) < 1.0e-10 && x[1]>=0.5;}) 
CELL_PREDICATE(BottomRightTest, {return fabs(x[0]-1.0) < 1.0e-10 && x[1]<=0.5;}) 

CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;}) 
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;}) 

CELL_PREDICATE(InterfaceTest, {return fabs(x[1]-0.5) < 1.0e-10;}) 

CELL_PREDICATE(BottomZoneTest, {return x[1] <= 0.5;}) 
CELL_PREDICATE(TopZoneTest, {return x[1] >= 0.5;}) 


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
      int nx = 32;
      int ny = 32;
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, npx, 
        0.0,  1.0, ny, npy, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter(); 
      CellFilter faces = new DimensionalCellFilter(1);
      
      CellFilter topZone = interior.subset(new TopZoneTest());
      CellFilter bottomZone = interior.subset(new BottomZoneTest());
      CellFilter interface = faces.subset(new InterfaceTest());
      CellFilter topLeft = faces.subset(new TopLeftTest());
      CellFilter bottomLeft = faces.subset(new BottomLeftTest());
      CellFilter topRight = faces.subset(new TopRightTest());
      CellFilter bottomRight = faces.subset(new BottomRightTest());

      CellFilter top = faces.subset(new TopPointTest());
      CellFilter bottom = faces.subset(new BottomPointTest());

      double mu1 = 0.5;
      double mu2 = 1.0;

      /* Create unknown function */
      Expr H1x = new UnknownFunction(new Lagrange(1), "H1x");
      Expr H1y = new UnknownFunction(new Lagrange(1), "H1y");
      Expr H1 = List(H1x, H1y);
      Expr B1 = mu1*H1;

      Expr H2x = new UnknownFunction(new Lagrange(1), "H2x");
      Expr H2y = new UnknownFunction(new Lagrange(1), "H2y");
      Expr H2 = List(H2x, H2y);
      Expr B2 = mu2*H2;
      
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr grad = List(dx, dy);
      
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      
      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(2);

      WatchFlag watch("watch me");
      watch.setParam("integration setup", 6);
      watch.deactivate();

      const double pi = 4.0*atan(1.0);
      double alpha1 = pi/4.0;
      double H1x0 = sin(alpha1);
      double H1y0 = cos(alpha1);

      double alpha2 = atan(mu2/mu1 * tan(alpha1));
      double absH2 = sqrt(pow(sin(alpha1),2) + pow(mu1/mu2,2)*pow(cos(alpha1),2));
      double H2x0 = absH2*sin(alpha2);
      double H2y0 = absH2*cos(alpha2);

      Out::os() << "exact soln (zone 1): {" << H1x0 << ", " << H1y0 << "}" << std::endl;
      Out::os() << "exact soln (zone 2): {" << H2x0 << ", " << H2y0 << "}" << std::endl;

      
      Expr sqResid = 
        Integral(bottomZone, curl(H1)*curl(H1) + div(B1)*div(B1), quad2)
        + Integral(topZone, curl(H2)*curl(H2) + div(B2)*div(B2), quad2)
        /* interface condition */
        + Integral(interface, (B2[1]-B1[1])*(B2[1]-B1[1]), quad2)
        + Integral(interface, (H2x-H1x)*(H2x-H1x), quad2, watch)
        /* BC on left and right surfaces */
        + Integral(topLeft+topRight, (dx*H2x)*(dx*H2x)+(dx*H2y)*(dx*H2y),
          quad2)
        + Integral(bottomLeft+bottomRight, (dx*H1x)*(dx*H1x)+(dx*H1y)*(dx*H1y),
          quad2)
        /* BC on top  */
        + Integral(top, (dy*H2x)*(dy*H2x)+(dy*H2y)*(dy*H2y),
          quad2)
        /* BC on bottom  */
        + Integral(bottom, (H1x-H1x0)*(H1x-H1x0) + (H1y-H1y0)*(H1y-H1y0),
          quad2);
      

      
      Functional sq(mesh, sqResid, vecType);
      
      Expr dum;
      Expr zero = new Sundance::ZeroExpr();
      Expr u = List(H1x, H1y, H2x, H2y);
      Expr v0 = List(zero,zero,zero,zero);
      LinearProblem prob = sq.linearVariationalProb(u, v0, u, dum, dum);

      ParameterXMLFileReader reader("aztec-ifpack.xml");

      ParameterList solverParams = reader.getParameters();

      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);


      Expr soln = prob.solve(solver);


      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("mag2D");
      w.addMesh(mesh);
      w.addField("H_1_x", new ExprFieldWrapper(soln[0]));
      w.addField("H_1_y", new ExprFieldWrapper(soln[1]));
      w.addField("H_2_x", new ExprFieldWrapper(soln[2]));
      w.addField("H_2_y", new ExprFieldWrapper(soln[3]));
      w.write();


      Out::os() << "exact soln (zone 1): {" << H1x0 << ", " << H1y0 << "}" << std::endl;
      Out::os() << "exact soln (zone 2): {" << H2x0 << ", " << H2y0 << "}" << std::endl;

      Expr errExpr_x_1 = Integral(bottom, pow(H1x0 - soln[0],2), quad2);
      Expr errExpr_x_2 = Integral(top, pow(H2x0 - soln[2],2), quad2);

      Expr errExpr_y_1 = Integral(bottom, pow(H1y0 - soln[1],2), quad2);
      Expr errExpr_y_2 = Integral(top, pow(H2y0 - soln[3],2), quad2);

      double errX1 = evaluateIntegral(mesh, errExpr_x_1);
      double errX2 = evaluateIntegral(mesh, errExpr_x_2);
      double errY1 = evaluateIntegral(mesh, errExpr_y_1);
      double errY2 = evaluateIntegral(mesh, errExpr_y_2);


      Out::os() << "zone 1, x error = " << sqrt(errX1) << std::endl;
      Out::os() << "zone 2, x error = " << sqrt(errX2) << std::endl;
      Out::os() << "zone 1, y error = " << sqrt(errY1) << std::endl;
      Out::os() << "zone 2, y error = " << sqrt(errY2) << std::endl;

      double errorSq = errX1 + errX2 + errY1 + errY2;
      double tol = 1.0e-6;
      Sundance::passFailTest(::sqrt(errorSq), tol);

    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}
  Sundance::finalize();
  return Sundance::testStatus(); 
}
