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

/** 
 * Solves the Poisson equation in 2D
 */


int main(int argc, char** argv)
{
  
  try
		{
      GlobalMPISession session(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedRectangleMesher. */
      MeshType meshType = new BasicSimplicialMeshType();

      MeshSource mesher 
        = new ExodusNetCDFMeshReader("../../../tests-std-framework/Problem/quarterCylinderCoarse.ncdf", meshType);
      Mesh mesh = mesher.getMesh();


      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter faces = new DimensionalCellFilter(2);
      CellPredicate topFunc = new LabelCellPredicate(1);
      CellPredicate bottomFunc = new LabelCellPredicate(2);
      CellPredicate xNormalFunc = new LabelCellPredicate(3);
      CellPredicate yNormalFunc = new LabelCellPredicate(5);
      CellFilter top = faces.subset(topFunc);
      CellFilter bottom = faces.subset(bottomFunc);
      CellFilter xNormalFace = faces.subset(xNormalFunc);
      CellFilter yNormalFace = faces.subset(yNormalFunc);

      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
			Expr vx = new TestFunction(new Lagrange(2), "vx");
			Expr vy = new TestFunction(new Lagrange(2), "vy");
			Expr vz = new TestFunction(new Lagrange(2), "vz");
			Expr ux = new UnknownFunction(new Lagrange(2), "ux");
			Expr uy = new UnknownFunction(new Lagrange(2), "uy");
			Expr uz = new UnknownFunction(new Lagrange(2), "uz");

      /* Create differential operator and coordinate functions */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr dz = new Derivative(2);
      Expr grad = List(dx, dy, dz);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      Expr z = new CoordExpr(2);

      
      /* Young's modulus */
      double E = 1.0;
      /* Poisson's ratio */
      double nu = 0.25;

      Expr lambda = E*nu/(1.0 + nu)/(1.0 - 2.0*nu);
      Expr mu = E/2.0/(1.0 + nu);
      Expr a = 2.0*mu + lambda;
      
      Expr D = List(List(a, lambda, lambda, 0.0, 0.0, 0.0),
                    List(lambda, a, lambda, 0.0, 0.0, 0.0),
                    List(lambda, a, lambda, 0.0, 0.0, 0.0),
                    List(0.0, 0.0, 0.0,     mu,  0.0, 0.0),
                    List(0.0, 0.0, 0.0,     0.0,  mu, 0.0),
                    List(0.0, 0.0, 0.0,     0.0, 0.0,  mu));

      Expr A = List(List(dx,  0.0, 0.0),
                    List(0.0,  dy, 0.0),
                    List(0.0, 0.0,  dz),
                    List(dy,   dx, 0.0),
                    List(dz,  0.0,  dx),
                    List(0.0,  dz,  dy));

      Expr strain = A*List(ux, uy, uz);
      Expr varStrain = A*List(vx, vy, vz);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(2);

      /* Define the weak form */
      Expr eqn = Integral(interior, varStrain*(D*strain), quad2)
        + Integral(top, -1.0*vz, quad2);

      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(xNormalFace, vx*ux, quad2)
        + EssentialBC(yNormalFace, vy*uy, quad2)
        + EssentialBC(bottom, vz*uz, quad2);


      Assembler::workSetSize() = 100;
      FunctionalEvaluator::workSetSize() = 100;

      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, List(vx, vy, vz), 
                         List(ux, uy, uz), vecType);

      /* Create an Aztec solver */
      std::map<int,int> azOptions;
      std::map<int,double> azParams;

      azOptions[AZ_solver] = AZ_gmres;
      azOptions[AZ_precond] = AZ_dom_decomp;
      azOptions[AZ_subdomain_solve] = AZ_ilu;
      azOptions[AZ_graph_fill] = 1;
      azParams[AZ_max_iter] = 1000;
      azParams[AZ_tol] = 1.0e-10;

      LinearSolver<double> solver = new AztecSolver(azOptions,azParams);

      Expr soln = prob.solve(solver);

      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("Elasticity3d");
      w.addMesh(mesh);
      w.addField("ux", new ExprFieldWrapper(soln[0]));
      w.addField("uy", new ExprFieldWrapper(soln[1]));
      w.addField("uz", new ExprFieldWrapper(soln[2]));
      w.write();


    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}
  TimeMonitor::summarize();
  
}
