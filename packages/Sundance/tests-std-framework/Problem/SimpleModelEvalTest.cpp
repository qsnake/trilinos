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

#include "SundanceDefs.hpp"
#include "SundanceDefs.hpp"



#ifdef HAVE_ENABLED_MOOCHO



#include "MoochoPack_MoochoThyraSolver.hpp"
#include "NLPInterfacePack_NLPFirstOrderThyraModelEvaluator.hpp"
#include "Thyra_DefaultModelEvaluatorWithSolveFactory.hpp"
#include "Thyra_DefaultRealLinearSolverBuilder.hpp"

#include "SundanceNLPModelEvaluator.hpp"
#include "Sundance.hpp"


CELL_PREDICATE_(LeftPointTest, {return fabs(x[0]) < 1.0e-10;});

CELL_PREDICATE_(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;});

  
namespace Thyra {


class SimpleSundanceModel : public SundanceNLPModelEvaluator
{
public:
  /** */
  SimpleSundanceModel(const VectorType<double>& vecType);
  
  /** */
  Array<double> exactSoln() const ;
  
  /** */
  Mesh mesh() const {return mesh_;}
  
private:
  Mesh mesh_;
};


} // namespace Thyra


int main(int argc, char** argv)
{
  using MoochoPack::MoochoSolver;
  using MoochoPack::MoochoThyraSolver;
	using NLPInterfacePack::NLPFirstOrderThyraModelEvaluator;

  try
  {

    Sundance::init(&argc, &argv);

    // We will do our linear algebra using Epetra
    VectorType<double> vecType = new EpetraVectorType();

    RCP<Thyra::SimpleSundanceModel> ssModel 
      = rcp(new Thyra::SimpleSundanceModel(vecType));
    RCP<Thyra::ModelEvaluator<double> > model = ssModel;

    // Get the LOWSF object from Stratimikos

    Thyra::DefaultRealLinearSolverBuilder lowsfBuilder;
    lowsfBuilder.setParameterList(Teuchos::parameterList());
    // ToDo: Pass in a real parameter list to select the linear solver and
    // preconditioner options!  This will just use the default which is
    // Amesos/Klu.

    RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory
      = createLinearSolveStrategy(lowsfBuilder);

    // Create the final Thyra::ModelEvaluator to be passed to MOOCHO

    RCP<Thyra::ModelEvaluator<double> > modelWithLOWS
      = rcp(
        new Thyra::DefaultModelEvaluatorWithSolveFactory<double>(
          model, lowsFactory
          )
        );

    // Create the solver object
    MoochoThyraSolver solver;
      
    // Set the Thyra model that defines the NLP
    solver.setModel(modelWithLOWS);
      
    // Solve the NLP
    const MoochoSolver::ESolutionStatus	solution_status = solver.solve();

    Array<double> exact = ssModel->exactSoln();
    Array<double> soln = ssModel->parameters();
    cout << "exact solution: " << exact << std::endl;
    cout << "numerical solution: " << soln << std::endl;
    double error = 0.0;
    for (int i=0; i<exact.size(); i++) error += pow(exact[i]-soln[i], 2.0);
    error = sqrt(error);

    bool OK = solution_status == MoochoSolver::SOLVE_RETURN_SOLVED;
    Sundance::passFailTest("MOOCHO converged", OK, error, 1.0e-3);

  }
	catch(std::exception& e)
  {
    Sundance::handleException(e);
  }

  Sundance::finalize(); return Sundance::testStatus(); 

}


namespace Thyra {


SimpleSundanceModel::SimpleSundanceModel(const VectorType<double>& vecType)
  : SundanceNLPModelEvaluator(vecType),
    mesh_()
{
  MeshType meshType = new BasicSimplicialMeshType();
  int np = MPIComm::world().getNProc();
  MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 10, meshType);
  mesh_ = mesher.getMesh();

  /* Create a cell filter that will identify the maximal cells
   * in the interior of the domain */
  CellFilter interior = new MaximalCellFilter();
  CellFilter points = new DimensionalCellFilter(0);
  CellFilter leftPoint = points.subset(new LeftPointTest());
  CellFilter rightPoint = points.subset(new RightPointTest());

  /* Create unknown and test functions, discretized using first-order
   * Lagrange interpolants */
  Expr u = new UnknownFunction(new Lagrange(2), "u");
  Expr v = new TestFunction(new Lagrange(2), "v");

  /* Create differential operator and coordinate function */
  Expr dx = new Sundance::Derivative(0);
  Expr x = new CoordExpr(0);

  const double pi = 4.0*atan(1.0);

  Expr source = 0.0;
  Array<Expr> a(4);

  for (int i=1; i<=a.size(); i++) 
  {
    a[i-1] = new Sundance::Parameter(1.0);
    source = source + pi*pi*i*i* a[i-1] * sin(i*pi*x);
  }
  Expr alpha = new ListExpr(a);

  /* We need a quadrature rule for doing the integrations */
  QuadratureFamily quad = new GaussianQuadrature(4);

      
  /* Define the weak form */
  Expr eqn = Integral(interior, (dx*v)*(dx*u) + v*source, quad);
  /* Define the Dirichlet BC */
  Expr bc = EssentialBC(leftPoint, v*u, quad)
    + EssentialBC(rightPoint, v*u, quad);

  /* Create a discrete space, and discretize the function 1.0 on it */
  DiscreteSpace discSpace(mesh_, new Lagrange(2), vecType);
  Expr u0 = new DiscreteFunction(discSpace, 1.0, "u0");

  /* create the forward problem */
  NonlinearProblem prob(mesh_, eqn, bc, v, u, u0, vecType);
      
  /* Write the sensitivity problems by hand. This will be unnecessary once
   * parametric differentiation is online. */
  Array<LinearProblem> sensProb(alpha.size());
  for (int i=0; i<alpha.size(); i++)
  { 
    double w = (i+1)*(i+1)*pi*pi;
    Expr sensEqn = Integral(interior, 
      (dx*v)*(dx*u) + w*v*sin((i+1)*pi*x), quad);
    Expr sensBC = EssentialBC(leftPoint, v*u, quad)
      + EssentialBC(rightPoint, v*u, quad);
    sensProb[i] = LinearProblem(mesh_, sensEqn, sensBC, v, u, vecType);
  }

  Expr sourceBasis = List(sin(pi*x), sin(2.0*pi*x), sin(3.0*pi*x), sin(4.0*pi*x));
  Expr uStar = 0.0;
  for (int i=0; i<exactSoln().size(); i++)
  {
    uStar = uStar - sourceBasis[i] * exactSoln()[i];
  }
    
  Expr objective = Integral(interior, 0.5*pow(u - uStar, 2.0), quad);
  Functional obj(mesh_, objective, vecType);

  initialize(alpha, u, u0, prob, sensProb, obj);
}


Array<double> SimpleSundanceModel::exactSoln() const
{
  return tuple(1.0, 0.3, 0.2, 0.1);
}


} // namespace Thyra



#else // HAVE_ENABLED_MOOCHO


#include <iostream>

int main()
{
  std::cout << "moocho not present: test INACTIVE" << std::endl;
}



#endif // HAVE_ENABLED_MOOCHO
