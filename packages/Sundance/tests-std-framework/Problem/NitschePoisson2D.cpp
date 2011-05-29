#include "Sundance.hpp"
#include "SundanceCellDiameterExpr.hpp"
#include "SundanceEvaluator.hpp"


NEW_CELL_PREDICATE(LeftPointTest)
{return fabs(x[0]) < 1.0e-10;}

NEW_CELL_PREDICATE(BottomPointTest)
{return fabs(x[1]) < 1.0e-10;}

NEW_CELL_PREDICATE(RightPointTest)
{return fabs(x[0]-1.0) < 1.0e-10;}

NEW_CELL_PREDICATE(TopPointTest)
{return fabs(x[1]-1.0) < 1.0e-10;}


void balanceXY(int n, int* npx, int* npy)
{
  int m = (int) floor(sqrt(n));
  for (int i=m; i>=1; i--)
  {
    if (n % i == 0) 
    {
      *npx = i;
      *npy = n/i;
      return ;
    }
  }

  *npx = n;
  *npy = 1;
}

/* weak form of poisson with Nitsche-type weak BC's */
Expr poissonEquationNitsche( bool splitBC, 
  Expr u ,
  Expr v ,
  Expr alpha ,
  QuadratureFamily quad )
{
  CellFilter interior = new MaximalCellFilter();
  CellFilter boundary = new BoundaryCellFilter();
  CellFilter left = boundary.subset( new LeftPointTest() );
  CellFilter right = boundary.subset( new RightPointTest() );
  CellFilter top = boundary.subset( new TopPointTest() );
  CellFilter bottom = boundary.subset( new BottomPointTest() );

  CellFilter allBdry = left+right+top+bottom;

  Expr dx = new Derivative(0);
  Expr dy = new Derivative(1);
  Expr x = new CoordExpr(0);
  Expr y = new CoordExpr(1);
  Expr grad = List( dx , dy );

  Expr uvTerm;
  if (splitBC)
  {
    Out::os() << "BC expressions split over domains" << std::endl;
    uvTerm = Integral( left , alpha*u * v , quad )
      + Integral( right , alpha*u * v , quad )
      + Integral( top , alpha*u * v , quad )
      + Integral( bottom , alpha*u * v , quad );
  }
  else
  {
    Out::os() << "BC expressions not split over domains" << std::endl;
    uvTerm = Integral( allBdry , alpha*u * v , quad );
  }
  

  const double pi = 4.0*atan(1.0);
  Expr force = 2.0*pi*pi*sin(pi*x)*sin(pi*y);
  return Integral( interior , (grad*v) * (grad*u) - force * v , quad )
    /* du/dn term */
    - Integral( left , -(dx*u)*v , quad )
    - Integral( top , (dy*u)*v , quad )
    - Integral( right , (dx*u)*v , quad )
    - Integral( bottom , -(dy*u)*v , quad )
    /* dv/dn term */
    - Integral( left , -(dx*v)*u , quad )
    - Integral( top , (dy*v)*u , quad )
    - Integral( right , (dx*v)*u , quad )
    - Integral( bottom , -(dy*v)*u , quad )
    /* u,v term  -- alpha = C / h */
    + uvTerm;
}


int main( int argc , char **argv )
{
  try {
    int nx = 128;
    double C = 4.0;
    std::string solverFile = "aztec-ml.xml";
    Sundance::setOption("nx", nx, "number of elements in x");
    Sundance::setOption("C", C, "Nitsche penalty");
    Sundance::setOption("solver", solverFile, "name of XML file for solver");
    
    Sundance::init( &argc , &argv );
    int np = MPIComm::world().getNProc();
    int npx = -1;
    int npy = -1;
    balanceXY(np, &npx, &npy);
    TEST_FOR_EXCEPT(npx < 1);
    TEST_FOR_EXCEPT(npy < 1);
    TEST_FOR_EXCEPT(npx * npy != np);

    VectorType<double> vecType = new EpetraVectorType();

    const int k = 1;
    const int splitBC = 1;

    MeshType meshType = new BasicSimplicialMeshType();
    MeshSource mesher = new PartitionedRectangleMesher( 0.0 , 1.0 , nx , npx ,
      0.0 , 1.0 , nx , npy,
      meshType );
    Mesh mesh = mesher.getMesh();

    BasisFamily L = new Lagrange( k );

    Expr u = new UnknownFunction( L , "u" );
    Expr v = new TestFunction( L , "v" );
    QuadratureFamily quad = new GaussianQuadrature( 2 * k );
    
    Expr h = new CellDiameterExpr();
    Expr alpha = C / h; 
    Expr eqn = poissonEquationNitsche( splitBC, u , v , alpha , quad );
    Expr bc;


    LinearProblem prob( mesh , eqn , bc , v , u , vecType);

#ifdef HAVE_CONFIG_H
    ParameterXMLFileReader reader(searchForFile("SolverParameters/" + solverFile));
#else
      ParameterXMLFileReader reader(solverFile);
#endif
    ParameterList solverParams = reader.getParameters();
    LinearSolver<double> solver 
      = LinearSolverBuilder::createSolver(solverParams);


    Expr soln = prob.solve( solver );


    FieldWriter w = new VTKWriter( "NitschePoisson2D" );
    w.addMesh( mesh );
    w.addField( "u" , new ExprFieldWrapper( soln ) );
    w.write();

    Expr x = new CoordExpr(0);
    Expr y = new CoordExpr(1);
    QuadratureFamily quad4 = new GaussianQuadrature(4);
    
    CellFilter interior = new MaximalCellFilter();
    const double pi = 4.0*atan(1.0);
    Expr exactSoln = sin(pi*x)*sin(pi*y);    
    Expr err = exactSoln - soln;
    Expr errExpr = Integral(interior, 
      err*err,
      quad4);
    FunctionalEvaluator errInt(mesh, errExpr);

    double errorSq = errInt.evaluate();
    cout << "error norm = " << sqrt(errorSq) << std::endl << std::endl;
    

    Sundance::passFailTest(sqrt(errorSq), 1.0e-4);
  }
  catch (std::exception &e) 
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); 
  return Sundance::testStatus();
}
