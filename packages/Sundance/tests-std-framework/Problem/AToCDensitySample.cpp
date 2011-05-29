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
#include "SundanceUserDefFunctor.hpp"
#include "SundanceUserDefOp.hpp"





/** 
 * 
 */

int main(int argc, char** argv)
{
  try
		{
      Sundance::init(&argc, &argv);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedLineMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(-1.0, 1.0, 32, 1,
                                                         -1.0, 1.0, 32, 1,
                                                         meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      BasisFamily L1 = new Lagrange(1);
      DiscreteSpace discSpace(mesh, List(L1, L1), vecType);
      
      /* Discretize some expression for the force. We'll pick a linear function
       * so that it can be interpolated exactly, letting us check the 
       * validity of our interpolations. */
      L2Projector proj(discSpace, List(x, y));
      Expr F = proj.project();

      
      /* create a sampler */
      cout << "making grid" << std::endl;
      AToCPointLocator locator(mesh, interior, createVector(tuple(200, 200)));
      
      AToCDensitySampler sampler(locator, vecType);

      CToAInterpolator forceInterpolator(locator, F);

      cout << "making points" << std::endl;
      /* create a bunch of particles */
      int nCells = mesh.numCells(2);
      int nPts = 15000;

      Array<double> pos(2*nPts);
      Array<double> f(F.size() * nPts);
      Array<Point> physPts;

      /* We'll generate random sample points in a way that lets us make an exact check
       * of the density recovery. We pick random cells, then random local coordinates
       * within each cell. This way, we can compute the density exactly as we
       * go, giving us something to check the recovered density against. */
      Array<int> counts(nCells);
      for (int i=0; i<nPts; i++)
        {
          /* pick a random cell */
          int cell = (int) floor(nCells * drand48());
          counts[cell]++;
          /* generate a point in local coordinates */
          double s = drand48();
          double t = drand48() * (1.0-s);
          Point refPt(s, t);
          /* map to physical coordinates */
          mesh.pushForward(2, tuple(cell), tuple(refPt), physPts);
          Point X = physPts[0];
          pos[2*i] = X[0];
          pos[2*i+1] = X[1];
        }

      cout << "sampling..." << std::endl;
      Expr density = sampler.sample(createVector(pos), 1.0);

      cout << "computing forces..." << std::endl;
      forceInterpolator.interpolate(pos, f);

      double maxForceErr = 0.0;
      for (int i=0; i<nPts; i++)
        {
          double x0 = pos[2*i];
          double y0 = pos[2*i+1];
          double fx = x0;
          double fy = y0;
          double df = ::fabs(fx - f[2*i]) + ::fabs(fy - f[2*i+1]);
          maxForceErr = max(maxForceErr, df);
        }
      cout << "max force error = " << maxForceErr << std::endl;

      cout << "writing..." << std::endl;

       /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("Density2d");
      w.addMesh(mesh);
      w.addField("rho", new ExprFieldWrapper(density));
      w.write();

      double errorSq = 0.0;

      double tol = 1.0e-6;
      Sundance::passFailTest(::sqrt(errorSq), tol);
    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}


