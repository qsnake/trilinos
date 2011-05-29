/* 
 * <Ignore> 
 * Copyright (2009) Kevin Long
 * Department of Mathematics and Statistics
 * Texas Tech University
 *
 * This library is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of the
 * License, or (at your option) any later version.
 *  
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *                                                                           
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA                                                                       
 * 
 *</Ignore>
 */

/* 
 * <Ignore> 
 * This boilerplate doesn't go in the documentation 
 */

#include "Sundance.hpp"
#include "SundanceElementIntegral.hpp"
using Sundance::List;

// This test depends on Exodus, so skip it if Expdus hasn't been enabled. 
#if defined(HAVE_SUNDANCE_EXODUS)

/* 
 * </Ignore> 
 */

/*
 * \documentclass[10pt]{report}
 * \usepackage{epsfig}
 * \begin{document}
 * 
 */


/* 
 * This example program sets up and solves the Poisson equation 
 * on
 * the domain shown in figure~\ref{F:PoissonDemoInputMesh}. The 
 *
 * \begin{figure}[h]
 * \begin{center}\includegraphics[width=4.0in]{plateWithHoleGeom.png}
 * \caption{
 * Initial mesh for 3D Poisson example.
 * }
 * \label{F:PoissonDemoInputMesh}
 * \end{center}\end{figure}
 */


/*
 *
 * <Header level="subsubsection" name="boilerplate">
 * C++ boilerplate
 * </Header>
 *
 * As always, a C++ program begins with \verb+main()+. We enclose the
 * body of the program inside a try/catch block for exception
 * handling. In most documentation, we'll ignore this boilerplate code.
 */

int main(int argc, char** argv)
{
  try
  {
    /*
     *
     * <Header level="subsubsection" name="init">
     * Initialization code
     * </Header>
     *
     * Every simulator starts with a call to \verb+Sundance::init()+, which
     * processes command-line arguments and initializes MPI. The
     * command-line handling is done using the Teuchos command line
     * parser. Code to set these options and their defaults
     * is normally the only work to be done before the call to \verb+init()+.
     */
    std::string meshFile="plateWithHole3D-1";
    std::string solverFile = "aztec-ml.xml";
    Sundance::setOption("meshFile", meshFile, "mesh file");
    Sundance::setOption("solver", solverFile, 
      "name of XML file for solver");

    /*
     * Now that the command-line parser has been set up, call \verb+init()+.
     */
    Sundance::init(&argc, &argv);

    // This next line is just a hack to deal with some 
    // transitional code in the
    // element integration logic. 
    Sundance::ElementIntegral::alwaysUseCofacets() = false;
    /* 
     * <Header level="subsubsection" name="vector_type">
     * Creation of vector type
     * </Header>
     * 
     * Next we create a \verb+VectorType+ object. A \verb+VectorType+
     * is an abstract factory which produces \verb+VectorSpace+ objects.
     * A \verb+VectorSpace+ is itself a factory which can produce 
     * \verb+Vector+ objects which are, logically enough, vectors.
     * 
     * Why this heirarchy of factories? Vectors may need to be created
     * inside code far below the user level, and the \verb+VectorSpace+ 
     * encapsulates the information needed to build a vector of a 
     * given type, size, and distribution over processors. The 
     * \verb+VectorType+ is needed because we might need to select a 
     * particular {\it software representation} for a vector, {\it e.g.}, 
     * Trilinos, Petsc, or some other library. 
     *
     * By using an \verb+EpetraVectorType+, we will be creating
     * \verb+EpetraVectorSpace+ vector spaces, which in turn
     * create \verb+EpetraVector+ vectors.
     */
    VectorType<double> vecType = new EpetraVectorType();

    /* 
     * <Header level="subsubsection" name="mesh">
     * Creation of mesh
     * </Header>
     * 
     * The creation of a mesh object involves several intermediate
     * objects: a \verb+MeshType+ and a \verb+MeshSource+. The 
     * \verb+MeshType+ specifies what mesh data structure will be used,
     * and the \verb+MeshSource+ builds that data structure. How it is
     * built will depend on the \verb+MeshSource+ subtype used. An
     * \verb+ExodusMeshReader+, for instance, reads from an Exodus II file,
     * while a \verb+PartitionedRectangleMesher+ builds a uniform 
     * triangulation of a rectangle.
     * 
     * Once the \verb+MeshType+ and \verb+MeshSource+ objects have been
     * defined, the mesh is created (read, built, whatever) by a call
     * to the \verb+getMesh()+ method of \verb+MeshSource+.
     */
    MeshType meshType = new BasicSimplicialMeshType();
      
    MeshSource meshSrc
      =  new ExodusMeshReader(meshFile, meshType);
    Mesh mesh = meshSrc.getMesh();

    /* 
     * <Header level="subsubsection" name="cell_filter">
     * Specification of geometric regions
     * </Header>
     * 
     * We'll need to specify subsets of the mesh on which equations
     * or boundary conditions are defined. In many FEA codes this is
     * done by explicit definition of element blocks, node sets, and
     * side sets. Rather than working with sets explicitly at the user
     * level, we instead work with filtering rules that produce 
     * sets of cells. These rules are represented by \verb+CellFilter+
     * objects. 
     *
     * \verb+MaximalCellFilter+
     * selects all cells having dimension equal to the spatial dimension
     * of the mesh. 
     */
    CellFilter interior = new MaximalCellFilter();

    
    /* 
     * \verb+DimensionalCellFilter+
     * selects all cells of a specified dimension. Here we select all 
     * edges. Boundary conditions will be applied on a subset of these.
     */
    CellFilter edges = new DimensionalCellFilter(2);

    CellFilter south = edges.labeledSubset(1);
    CellFilter east = edges.labeledSubset(2);
    CellFilter north = edges.labeledSubset(3);
    CellFilter west = edges.labeledSubset(4);
    CellFilter hole = edges.labeledSubset(5);
    CellFilter down = edges.labeledSubset(6);
    CellFilter up = edges.labeledSubset(7);

    /* 
     * <Header level="subsubsection" name="symb_setup">
     * Setup of symbolic problem description
     * </Header>
     * 
     * Create unknown and test functions discretized on the space
     * first-order Lagrange polynomials. 
     */
    BasisFamily basis = new Lagrange(1);
    Expr u = new UnknownFunction(basis, "u");
    Expr v = new TestFunction(basis, "v");

    /* 
     * Create differential operators and coordinate functions. Directions
     * are indexed starting from zero. The \verb+List()+ function can 
     * collect expressions into a vector. 
     */
    Expr dx = new Derivative(0);
    Expr dy = new Derivative(1);
    Expr dz = new Derivative(2);
    Expr grad = List(dx, dy, dz);

    Expr x = new CoordExpr(0);
    Expr y = new CoordExpr(1);
    Expr z = new CoordExpr(2);

    /* 
     * We need a quadrature rule for doing the integrations 
     */
    QuadratureFamily quad2 = new GaussianQuadrature(2);
    QuadratureFamily quad4 = new GaussianQuadrature(4);

    /* 
     * Having defined the necessary symbolic expressions, cell filters, and
     * quadrature rules, we can now use them to create the weak form  
     * plus essential BCs. 
     *
     * We also introduce one additional feature that aids debugging: 
     * a watch flag. Any term (or collection of terms) in a weak form
     * or essential BC can be marked with a watch flag which can provide
     * fine-tuned control over the diagnostic output to be produced. 
     * Generally when debugging, you will have isolated the problem to
     * one or two terms in an expression, and the watch flag lets you 
     * eliminate the distracting output from other terms. 
     * 
     * A watch flag can be disabled by invoking the \verb+deactivate()+
     * method. Verbosity levels for various stages of the calculation
     * can be set by the \verb+setParam()+ method; increasing the
     * verbosity level increases the amount of information printed.
     */
    WatchFlag watchMe("watch me");
    watchMe.setParam("symbolic preprocessing", 1);
    watchMe.setParam("discrete function evaluation", 3);
    watchMe.setParam("integration setup", 6);
    watchMe.setParam("integration", 6);
    watchMe.setParam("fill", 6);
    watchMe.setParam("evaluation", 2);
    watchMe.deactivate();  // Deactive
    Expr source = 0.0;
    Expr eqn 
      = Integral(interior, (grad*u)*(grad*v), quad2, watchMe);

    /* 
     * <Header level="subsubsection" name="essential_bc">
     * Setup of essential boundary conditions
     * </Header>
     */
    Expr h = new CellDiameterExpr();
    Expr bc = EssentialBC(west, v*u/h, quad2) +
      EssentialBC(east, v*(u-1.0)/h, quad2);

    /* 
     * <Header level="subsubsection" name="lin_prob">
     * Creation of linear problem
     * </Header>
     *
     * We have now defined the mesh, equations, boundary conditions, and
     * are ready to produce a discrete problem. 
     *
     * Note: \verb+LinearProblem+ is a lightweight user-level interface
     * to an \verb+Assembler+ object which actually does the work of
     * building matrices and vectors. \verb+Assembler+ is also used
     * under the hood for the assembly of Jacobians and residuals for
     * nonlinear problems and for the calculation of functional
     * values and gradients. \verb+LinearProblem+ ensures that the
     * \verb+Assembler+ is constructed properly, controls the call to
     * \verb+Assembler+ for building the matrix and vector, invokes 
     * the linear solver and checks for convergence, and wraps the 
     * solution vector in a \verb+DiscreteFunction+ object so that it
     * can be used in symbolic specification of future problems.
     *
     * {\bf TODO:} As I was writing this, I realized that LP is the logical
     * location for sanity checks on equation sets. I'd tried to
     * build that into \verb+EquationSet+, but it was complicated by
     * the need to deal with many problem types (with/without test functions,
     * and so on). Look into this.
     *
     * The simplest \verb+LinearProblem+ constructor accepts a 
     * mesh, the weak equations and essential boundary conditions, 
     * the test and unknown functions, and the vector type to be used.
     */
    LinearProblem prob(mesh, eqn, bc, v, u, vecType);

    /*
     * Other \verb+LinearProblem+ constructors will be 
     * required for more complicated cases:
     * \begin{itemize}
     * \item Forming block matrix systems ({\it e.g.}, segregated
     * preconditioners such as Kay-Lohgin or Elman-Howle-Tuminaro 
     * approximate Schur complement methods for saddle-point problems).
     * \item Holding a subset of variables fixed during one phase
     * of a calculation. This is rarely done at the user level; it is
     * usually done internally when setting up variations of a functional
     * in a reduced-space optimization problem.
     * \end{itemize}
     * These will be shown in subsequent examples.
     * 
     * <Header level="subsubsection" name="lin_solve">
     * Solution of linear problem
     * </Header>
     * 
     */
    
    ParameterXMLFileReader reader(solverFile);
    ParameterList solverParams = reader.getParameters();
    LinearSolver<double> solver 
      = LinearSolverBuilder::createSolver(solverParams);

    Expr soln = prob.solve(solver);

    /* 
     * <Header level="subsubsection" name="vizout">
     * Visualization output
     * </Header>
     * Output to storage or visualization files is done using an 
     * abstract \verb+FieldWriter+ interface. \verb+FieldWriter+
     * concrete types include \verb+VTKWriter+ for output to 
     * unstructured VTK format,  \verb+MatlabWriter+ for output
     * of simple column-oriented files for 1D problems, and
     * \verb+VerboseFieldWriter+ for a verbose dump showing facet
     * and cofacet relations.
     *
     * In this example, we use \verb+VTKWriter+. The solution as rendered
     * by Paraview is shown in figure~\ref{F:PoissonDemoSoln}.
     * \begin{figure}[h]
     * \begin{center}\includegraphics[width=4.0in]{plateWithHole3D.jpg}
     * \caption{
     * Solution to 3D Poisson example problem.
     * }
     * \label{F:PoissonDemoSoln}
     * \end{center}\end{figure}
     * 
     */
    FieldWriter w = new VTKWriter("PoissonDemo3D");
    w.addMesh(mesh);
    w.addField("soln", new ExprFieldWrapper(soln[0]));
    w.write();


    /* 
     * <Header level="subsubsection" name="postproc">
     * Postprocessing
     * </Header>
     *
     * Postprocessing can be done using the same symbolic language
     * as was used for the problem specification. Here, we define
     * an integral giving the flux, then evaluate it on the mesh. 
     */
    Expr n = CellNormalExpr(3, "n");
    Expr fluxExpr 
      = Integral(east + west, (n*grad)*soln, quad2); 
    double flux = evaluateIntegral(mesh, fluxExpr);
    Out::os() << "numerical flux = " << flux << std::endl;

    /* 
     * Let's compute a few other quantities, such as the centroid of
     * the mesh:
     */
    Expr xCMExpr = Integral(interior, x, quad2);
    Expr yCMExpr = Integral(interior, y, quad2);
    Expr zCMExpr = Integral(interior, z, quad2);
    double xCM = evaluateIntegral(mesh, xCMExpr);
    double yCM = evaluateIntegral(mesh, yCMExpr);
    double zCM = evaluateIntegral(mesh, zCMExpr);
    Out::os() << "centroid = (" << xCM << ", " << yCM 
              << ", " << zCM << ")" << std::endl;

    /*
     * Next, compute the first Fourier sine 
     * coefficient of the solution on the
     * hole.
     */
    Expr r = sqrt(x*x + y*y);
    Expr sinPhi = y/r;
    Expr fourierSin1Expr = Integral(hole, sinPhi*soln, quad2);
    double fourierSin1 = evaluateIntegral(mesh, fourierSin1Expr);
    Out::os() << "fourier sin m=1 = " << fourierSin1 << std::endl;

    /*
     * Check that the flux is acceptably close to zero. This
     * is just a sanity check to ensure the code doesn't get completely 
     * broken after a change to the library. 
     */
    Sundance::passFailTest(fabs(flux), 1.0e-3);

    /*
     * <Header level="subsubsection" name="finalize">
     * Finalization boilerplate
     * </Header>
     * Finally, we have boilerplate code for exception handling
     * and finalization. 
     */

  }
	catch(std::exception& e) 
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); return Sundance::testStatus(); 

  return Sundance::testStatus();
}

/*
 * All done!
 */

/* <Ignore> */
#else


int main(int argc, char** argv)
{
  Sundance::init(&argc, &argv);
  std::cout << "dummy PoissonDemo3D PASSED. Enable exodus to run the actual test" << std::endl;
  Sundance::finalize();
  return 0;
}


#endif

/* </Ignore> */

/* \end{document} */
