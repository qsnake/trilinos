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
#include "SundanceCFMeshPair.hpp"
#include "SundanceDOFMapBuilder.hpp"
#include "SundanceInhomogeneousNodalDOFMap.hpp"


/** 
 * Tests logical operations on cell filters
 */


CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;});

CELL_PREDICATE(ATest, {return x[0] <= 0.4;});

CELL_PREDICATE(BTest, {return x[0] >= 0.4 && x[0] <= 0.6;});

CELL_PREDICATE(CTest, {return x[0] >= 0.6;});
  
int main(int argc, void** argv)
{
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedLineMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      int nx = 10;
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx, meshType);
      Mesh mesh = mesher.getMesh();

      Expr x = new CoordExpr(0);

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter points = new DimensionalCellFilter(0);
      CellFilter leftPoint = points.subset(new LeftPointTest());
      CellFilter A = interior.subset(new ATest());
      CellFilter B = interior.subset(new BTest());
      CellFilter C = interior.subset(new CTest());

      Expr u1 = new UnknownFunction(new Lagrange(1));
      Expr u2 = new UnknownFunction(new Lagrange(1));

      Expr v1 = new TestFunction(new Lagrange(1));
      Expr v2 = new TestFunction(new Lagrange(1));

      QuadratureFamily quad = new GaussianQuadrature(2);
      Expr eqn = Integral(interior, v1*(u1 - 2.0), quad) 
        + Integral(A, v2*(u2 - x*u1), quad) 
        + Integral(B, v2*(u2 - 0.4*u1), quad) ;
      Expr bc;

      LinearProblem prob(mesh, eqn, bc, List(v1, v2), List(u1, u2), vecType);


      cout << "-------- ROW MAP --------------------------" << std::endl;
      prob.rowMap(0)->print(cout);

      cout << "-------- COL MAP --------------------------" << std::endl;
      prob.colMap(0)->print(cout);

#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/aztec.xml"));
#else
      ParameterXMLFileReader reader("aztec.xml");
#endif
      ParameterList solverParams = reader.getParameters();
      cout << "params = " << solverParams << std::endl;


      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      cout << "matrix = " << std::endl << prob.getOperator() << std::endl;
      cout << "RHS = " << std::endl << prob.getRHS() << std::endl;

      Expr soln = prob.solve(solver);

      Vector<double> vec = DiscreteFunction::discFunc(soln)->getVector();

      cout << "solution = " << vec << std::endl;
      
    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
