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



/**
 * This program tests the ability to compute derivatives of discrete
 * functions in thre dimensions.
 */

int main(int argc, char** argv)
{
  try
  {
   
    std::string meshFile="plate3D-0";
    int order = 2;
    Sundance::setOption("meshFile", meshFile, "mesh file");
    Sundance::setOption("order", order, "order");

    Sundance::init(&argc, &argv);

    // This next line is just a hack to deal with some 
    // transitional code in the
    // element integration logic. 
    Sundance::ElementIntegral::alwaysUseCofacets() = false;


    VectorType<double> vecType = new EpetraVectorType();

    MeshType meshType = new BasicSimplicialMeshType();
      
    MeshSource meshSrc
      =  new ExodusMeshReader(meshFile, meshType);
    Mesh mesh = meshSrc.getMesh();

    CellFilter interior = new MaximalCellFilter();

    
    BasisFamily basis = new Lagrange(order);
    DiscreteSpace discSpace(mesh, basis, vecType);


    Expr dx = new Derivative(0);
    Expr dy = new Derivative(1);
    Expr dz = new Derivative(2);
    Expr grad = List(dx, dy, dz);

    Expr x = new CoordExpr(0);
    Expr y = new CoordExpr(1);
    Expr z = new CoordExpr(2);
    

    /* */
    Expr f = 0.5*(x*x + y*y + z*z);
    L2Projector proj(discSpace, f);
    Expr u0 = proj.project();
    L2Projector proj2(discSpace, dx*u0 - dx*f);
    Expr du0_dx = proj2.project();
    L2Projector proj3(discSpace, dy*u0 - dy*f);
    Expr du0_dy = proj3.project();
    L2Projector proj4(discSpace, dz*u0 - dz*f);
    Expr du0_dz = proj4.project();

    Expr F = List(x,y,z);
    Expr n = CellNormalExpr(3, "n");
    Expr v = new TestFunction(basis);
    Expr u = new UnknownFunction(basis);
    CellFilter bdry = new BoundaryCellFilter();
    QuadratureFamily quad = new GaussianQuadrature(4);
    
    Expr eqn = Integral(interior, u*v + F*(grad*v) + v*(grad*F), quad)
      - Integral(bdry, v*(F*n), quad);
    Expr bc;
    LinearProblem prob(mesh, eqn, bc, v, u, vecType);
    
    ParameterXMLFileReader reader("amesos.xml");
    ParameterList solverParams = reader.getParameters();
    LinearSolver<double> solver 
      = LinearSolverBuilder::createSolver(solverParams);

    Expr soln = prob.solve(solver);

    FieldWriter w = new VTKWriter("DiscFunc3D");
    w.addMesh(mesh);
    w.addField("u", new ExprFieldWrapper(u0));
    w.addField("soln", new ExprFieldWrapper(soln));
    w.addField("u_x", new ExprFieldWrapper(du0_dx));
    w.addField("u_y", new ExprFieldWrapper(du0_dy));
    w.addField("u_z", new ExprFieldWrapper(du0_dz));
    w.write();

    Expr e0 = Integral(interior, soln, quad);
    Expr e1 = Integral(interior, du0_dx, quad);
    Expr e2 = Integral(interior, du0_dy, quad);
    Expr e3 = Integral(interior, du0_dz, quad);
    double err0 = fabs(evaluateIntegral(mesh, e0));
    double err1 = fabs(evaluateIntegral(mesh, e1));
    double err2 = fabs(evaluateIntegral(mesh, e2));
    double err3 = fabs(evaluateIntegral(mesh, e3));

    Out::os() << "error in Green's formula = " << err0 << std::endl;
    Out::os() << "error in du/dx = " << err1 << std::endl;
    Out::os() << "error in du/dy = " << err2 << std::endl;
    Out::os() << "error in du/dz = " << err3 << std::endl;

    
    Sundance::passFailTest(err0 + err1 + err2 + err3, 1.0e-10);

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
  std::cout << "dummy DiscFunc3D PASSED. Enable exodus to run the actual test" << std::endl;
  Sundance::finalize();
  return 0;
}


#endif

/* </Ignore> */

/* \end{document} */
