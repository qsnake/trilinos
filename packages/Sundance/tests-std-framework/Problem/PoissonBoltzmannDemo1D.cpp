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


CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})

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
     * Initialization code
     */
    std::string solverFile = "nox-amesos.xml";
    Sundance::setOption("solver", solverFile, 
      "name of XML file for solver");

    Sundance::init(&argc, &argv);

    // This next line is just a hack to deal with some 
    // transitional code in the
    // element integration logic. 
    Sundance::ElementIntegral::alwaysUseCofacets() = false;

    /* 
     * Creation of vector type
     */
    VectorType<double> vecType = new EpetraVectorType();

    /* 
     * Creation of mesh
     */
    MeshType meshType = new BasicSimplicialMeshType();
      
    MeshSource meshSrc
      = new PartitionedLineMesher(0.0, 1.0, 1000, meshType);
    Mesh mesh = meshSrc.getMesh();

    /* 
     * Specification of cell filters
     */
    CellFilter interior = new MaximalCellFilter();
    CellFilter pts = new DimensionalCellFilter(0);
    CellFilter left = pts.subset(new LeftPointTest());
    CellFilter right = pts.subset(new RightPointTest());
    

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
    Expr grad = List(dx);

    Expr x = new CoordExpr(0);

    /* 
     * We need a quadrature rule for doing the integrations 
     */
    QuadratureFamily quad2 = new GaussianQuadrature(2);
    QuadratureFamily quad4 = new GaussianQuadrature(4);

    /* 
     * Create the weak form and the BCs
     */
    Expr source = exp(-u/10.0);
    Expr eqn 
      = Integral(interior, (grad*u)*(grad*v)+source*v, quad2);

    Expr bc = EssentialBC(left, v*u, quad2) +
      EssentialBC(right, v*(u-1.0), quad2);

    /* 
     * <Header level="subsubsection" name="lin_prob">
     * Creation of initial guess
     * </Header>
     *
     * So far the setup has been almost identical to that for the linear
     * problem, the only difference being the nonlinear term in the
     * equation set. 
     */
    DiscreteSpace discSpace(mesh, basis, vecType);
    L2Projector proj(discSpace, 0.0);
    Expr u0 = proj.project();

    /* 
     * <Header level="subsubsection" name="lin_prob">
     * Creation of nonlinear problem
     * </Header>
     *
     * Similar to the setup of a \verb+LinearProblem+, the equation, BCs,
     * and mesh are put into a \verb+NonlinearProblem+ object which
     * controls the construction of the \verb+Assembler+ and its use
     * in building Jacobians and residuals during a nonlinear solve.
     */
    NonlinearProblem prob(mesh, eqn, bc, v, u, u0, vecType);

    /*
     *
     */
    
    ParameterXMLFileReader reader(solverFile);
    ParameterList solverParams = reader.getParameters();
    NOXSolver solver(solverParams); 
    
    prob.solve(solver);
    
    /* 
     * Visualization output
     */
    FieldWriter w = new MatlabWriter("PoissonBoltzmannDemo1D.dat");
    w.addMesh(mesh);
    w.addField("soln", new ExprFieldWrapper(u0));
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
    Expr n = CellNormalExpr(1, "n");
    Expr fluxExpr 
      = Integral(left+right, (n*grad)*u0, quad2); 
    double flux = evaluateIntegral(mesh, fluxExpr);
    Out::os() << "numerical flux = " << flux << std::endl;
    Expr sourceExpr 
      = Integral(interior, exp(-u0/10.0), quad2); 
    double src = evaluateIntegral(mesh, sourceExpr);
    Out::os() << "numerical src = " << src << std::endl;



    /*
     * Check that the flux is acceptably close to zero. This
     * is just a sanity check to ensure the code doesn't get completely 
     * broken after a change to the library. 
     */
    Sundance::passFailTest(fabs(flux-src), 1.0e-3);

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
