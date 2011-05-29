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

class BesselJFunc : public PointwiseUserDefFunctor0
{
public: 
  /** */
  BesselJFunc(int n) : PointwiseUserDefFunctor0("J_" + Teuchos::toString(n), 1, 1), 
                       n_(n) {;}

  /** */
  void eval0(const double* vars, double* f) const ;

private:
  int n_;
};

Expr BesselJ(int n, const Expr& r)
{
  return  new UserDefOp(r, rcp(new BesselJFunc(n)));
}


/** 
 * To plug in a user-defined function, you write a functor object,
 * deriving from UserDefFunctor, and implement your function in 
 * the eval() method of the functor. Any constant parameters (i.e., 
 * NOT expr objects) to your function, you can pass in as ctor 
 * arguments. 
 *
 * A user-defined functor will be accessed through a UserDefOp
 * object, whose ctor takes Sundance expression specifying the 
 * arguments of the function. This expression can be list-valued
 * in the event that your function is of more than one variable.
 *
 * The eval() method evaluates the function at a single point
 * in space. The argument to eval() is an array of doubles, the elements
 * of which are the values of the user-level argument Expr at
 * the current evaluation point. 
 */
class DrumFuncFunctor : public PointwiseUserDefFunctor0
{
public:
  /** */
  DrumFuncFunctor(const double& A, int n) 
    : PointwiseUserDefFunctor0("Drum_" + Teuchos::toString(n), 2, 1), A_(A), n_(n) {;}

  /** */
  void eval0(const double* vars, double* f) const ;

private:
  double A_;
  int n_;
};


/* Definition of eval() function. This is where the "meat" is. */

void DrumFuncFunctor::eval0(const double* vars, double* f) const
{
  double x = vars[0];
  double y = vars[1];

  double theta = atan2(y,x);
  double r = sqrt(x*x + y*y);

  f[0] = A_*cos(n_*theta)*jn(n_, r);
}



/* A little wrapper to be able to build your user-defined function
* in pretty syntax. */
Expr DrumFunc(int n, double A, const Expr& x, const Expr& y)
{
  return new UserDefOp(List(x,y), rcp(new DrumFuncFunctor(A, n)));
}



#define LOUD()                                          \
  {                                                     \
    verbosity<Evaluator>() = 5;               \
    verbosity<SparsitySuperset>() = 0;         \
    verbosity<EvalVector>() = 0;               \
    verbosity<EvaluatableExpr>() = 5;         \
    verbosity<AbstractEvalMediator>() = 5;    \
  }

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
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 10.0, 64, 1,
                                                         0.0, 10.0, 64, 1,
                                                         meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();


      DiscreteSpace discSpace(mesh, new Lagrange(2), vecType);

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      Expr r = sqrt(x*x + y*y);

      int n = 2;
      Expr Jn = BesselJ(n, r);

      Expr Dn = DrumFunc(n, 1.0, x, y);

      Expr JnDisc = L2Projector(discSpace, Jn).project();
      Expr DnDisc = L2Projector(discSpace, Dn).project();

      Expr errExpr = Integral(interior, 
                              0.01 * pow(JnDisc*(x*x - y*y)/r/r -DnDisc, 2.0),
                              new GaussianQuadrature(4) );

      double errorSq = evaluateIntegral(mesh, errExpr);

       /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("Bessel2d");
      w.addMesh(mesh);
      w.addField("J_n", new ExprFieldWrapper(JnDisc));
      w.addField("Drum_n", new ExprFieldWrapper(DnDisc));
      w.write();

      double tol = 4.0e-6;
      Sundance::passFailTest(::sqrt(errorSq), tol);
    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}



/* Definition of BesselJFunc::eval() */

void BesselJFunc::eval0(const double* vars, double* f) const
{
  f[0] = jn(n_, vars[0]);
}


