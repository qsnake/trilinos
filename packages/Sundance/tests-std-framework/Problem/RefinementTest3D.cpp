#include "Sundance.hpp"
#include "SundanceRivaraMesh.hpp"
#include "SundanceRivaraDriver.hpp"
#include "SundanceBubble.hpp"
#include "SundanceExodusWriter.hpp"

using namespace Sundance::Rivara;

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]+1.0) < 1.0e-10;}) 
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})

Expr poissonSolve(
  const Mesh& mesh, 
  const VectorType<double>& vecType,
  const LinearSolver<double>& solver
  );

Expr errEstimate(
  const Mesh& mesh, 
  const VectorType<double>& vecType,
  const Expr& soln,
  double& errNorm
  );

double scaleLength() {return 0.2;}

Expr density(const Expr& x)
{
  double a2 = std::pow(scaleLength(), 2.0);
  Expr x2 = x*x;
  Expr d = x2 + a2;
  return a2*(8.0*x2/d - 2.0)/d/d;
}

Expr exactSoln(const Expr& x)
{
  double a2 = std::pow(scaleLength(), 2.0);
  Expr x2 = x*x;
  Expr d = x2 + a2;
  return -a2/d + a2/(a2 + 1.0);
}

int main(int argc, char** argv)
{

  try
		{
			Sundance::init(&argc, &argv);

      /* Get a mesh */
      int nx = 32;
      int ny = 32;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(
        -1.0, 1.0, nx, 1, 
        -1.0,  1.0, ny, 1, meshType);

      Mesh mesh = mesher.getMesh();

      MeshTransformation extruder 
        = new ExtrusionMeshTransformation(0.0, 0.08, 1, meshType);

      mesh = extruder.apply(mesh);

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      VectorType<double> vecType = new EpetraVectorType();



#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/amesos.xml"));
#else
      ParameterXMLFileReader reader("amesos.xml");
#endif
      ParameterList solverParams = reader.getParameters();
      LinearSolver<double> solver
        = LinearSolverBuilder::createSolver( solverParams );

      double errEst = 1.0;
      double localTol = 2.0e-5;
      double errNorm = 1.0;

      CellFilter interior = new MaximalCellFilter();
      
      QuadratureFamily quad4 = new GaussianQuadrature(4);
      
      for (int i=0; i<10; i++)
      {
        Out::os() << "mesh has " << mesh.numCells(0) << " nodes " << std::endl;
        Out::os() << "solving" << std::endl;
        Expr soln = poissonSolve(mesh, vecType, solver);         
        Out::os() << "computing error estimate" << std::endl;
        Expr err = errEstimate(mesh, vecType, soln, errEst);
        DiscreteSpace ds(mesh, new Lagrange(0), vecType);
        L2Projector pr(ds, soln - exactSoln(new CoordExpr(0)));
        Expr exactDiff = pr.project();

        Out::os() << "error estimate = " << errEst << std::endl;
        errNorm = evaluateIntegral(mesh, Integral(interior, exactDiff*exactDiff, quad4));
        errNorm = std::sqrt(errNorm);
        Out::os() << "error norm = " << errNorm << std::endl;

        FieldWriter w = new VTKWriter("refined-" + Teuchos::toString(i));
        w.addMesh(mesh);
        w.addField("soln", new ExprFieldWrapper(soln));
        w.addField("err est", new ExprFieldWrapper(err));
        w.addField("err ex", new ExprFieldWrapper(exactDiff));
        w.write();

        Out::os() << "refining" << std::endl;
        RefinementTransformation ref(meshType, err, localTol, 1.0e-12);
        mesh = ref.apply(mesh);
        int numRefined = ref.numRefined();
        if (numRefined==0) break;
      }


      Sundance::passFailTest(errNorm, 2.0e-3);
    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}
	Sundance::finalize();
  return Sundance::testStatus(); 
}





Expr poissonSolve(
  const Mesh& mesh, 
  const VectorType<double>& vecType,
  const LinearSolver<double>& solver
  )
{
  Expr x = new CoordExpr(0);

  Expr grad = gradient(mesh.spatialDim());

  Expr u = new UnknownFunction(new Lagrange(1));
  Expr v = new TestFunction(new Lagrange(1));

  CellFilter interior = new MaximalCellFilter();
  CellFilter edges = new DimensionalCellFilter(mesh.spatialDim()-1);
  CellFilter left = edges.subset(new LeftPointTest());
  CellFilter right = edges.subset(new RightPointTest());

  QuadratureFamily quad4 = new GaussianQuadrature(4);

  Expr eqn = Integral(interior, (grad*u)*(grad*v) - v*density(x), quad4);

  Expr bc = EssentialBC(left+right, v*u,quad4);

  LinearProblem prob(mesh, eqn, bc, v, u, vecType);

  return prob.solve(solver);

}

Expr errEstimate(
  const Mesh& mesh, 
  const VectorType<double>& vecType,
  const Expr& u0,
  double& errNorm
  )
{
  Expr x = new CoordExpr(0);

  Expr grad = gradient(mesh.spatialDim());

  Expr v = new TestFunction(new Bubble(3));
  Expr u = new UnknownFunction(new Lagrange(0));

  CellFilter interior = new MaximalCellFilter();

  QuadratureFamily quad4 = new GaussianQuadrature(4);

  Expr residEqn 
    = Integral(interior, u*v+(grad*u0)*(grad*v) - v*density(x), quad4);
  Expr normEqn = Integral(interior, u*v + v, quad4);

  Expr bc;

  LinearProblem prob(mesh, residEqn, bc, v, u, vecType);
  Vector<double> eVec = prob.getRHS().abs();

  Expr err = prob.formSolutionExpr(eVec);

  errNorm = evaluateIntegral(mesh, Integral(interior, err[0], quad4));

  return err;
}

