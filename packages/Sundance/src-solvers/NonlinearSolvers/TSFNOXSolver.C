// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "TSFNOXSolver.H"         
#include "NOX_StatusTest_SafeCombo.H"         
#include "NOX.H"         
//#include "NOX_Parameter_Teuchos2NOX.H"         
#include "TSFLinearSolverBuilder.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceExceptions.hpp"
#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorImpl.hpp"
#include "TSFLinearOperatorImpl.hpp"
#endif


using namespace NOX;
using namespace NOX::TSF;
using namespace Teuchos;
using namespace TSFExtended;


static Time& noxSolverTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("NOX solve"); 
  return *rtn;
}

NOXSolver::NOXSolver(const ParameterList& params)
  : linSolver_(),
    statusTest_(),
    params_(),
    printParams_()
{
  TEST_FOR_EXCEPTION(!params.isSublist("NOX Solver"), runtime_error,
                     "did not find NOX Solver sublist in " << params);
  
  params_ = params.sublist("NOX Solver");
  /* NOX wants to have the process ID in a parameter list???? */
  params_.sublist("Printing").set("MyPID", MPIComm::world().getRank());

  if (params_.isSublist("Status Test"))
    {
      statusTest_ = StatusTestBuilder::makeStatusTest(params_);
    }
  else
    {
      RCP<StatusTest::Generic> A = rcp(new StatusTest::NormF(1.0e-12));
      RCP<StatusTest::Generic> B = rcp(new StatusTest::MaxIters(20));
      statusTest_ = 
        rcp(new StatusTest::SafeCombo(StatusTest::SafeCombo::OR, A, B));
    }
  
  if (params_.isSublist("Linear Solver"))
    {
      linSolver_ = LinearSolverBuilder::createSolver(params_);
    }
  else
  {
    TEST_FOR_EXCEPTION(!params_.isSublist("Linear Solver"),
      RuntimeError, "no linear solver specified in NOX parameters");
  }
  
  if (params_.isSublist("Printing"))
    {
      printParams_ = params_.sublist("Printing");
    }
  
  TEST_FOR_EXCEPTION(linSolver_.ptr().get()==0, runtime_error,
                     "null linear solver object in NOXSolver ctor");

  TEST_FOR_EXCEPTION(statusTest_.get()==0, runtime_error,
                     "null status test object in NOXSolver ctor");

}

NOXSolver::NOXSolver(const ParameterList& nonlinParams,
      const LinearSolver<double>& linSolver)
  : linSolver_(linSolver),
    statusTest_(),
    params_(),
    printParams_()
{
  Tabs tab(0);
  TEST_FOR_EXCEPTION(!nonlinParams.isSublist("NOX Solver"), runtime_error,
                     "did not find NOX Solver sublist in " << nonlinParams);
  
  params_ = nonlinParams.sublist("NOX Solver");
  /* NOX wants to have the process ID in a parameter list???? */
  params_.sublist("Printing").set("MyPID", MPIComm::world().getRank());

  if (params_.isSublist("Status Test"))
    {
      statusTest_ = StatusTestBuilder::makeStatusTest(params_);
    }
  else
    {
      RCP<StatusTest::Generic> A = rcp(new StatusTest::NormF(1.0e-12));
      RCP<StatusTest::Generic> B = rcp(new StatusTest::MaxIters(20));
      statusTest_ = 
        rcp(new StatusTest::SafeCombo(StatusTest::SafeCombo::OR, A, B));
    }
  
  if (params_.isSublist("Linear Solver"))
    {
      Out::root() << tab << "WARNING: linear solver in NOX parameter list "
        "will be overridden by alternate solver" << std::endl;
    }
  
  if (params_.isSublist("Printing"))
    {
      printParams_ = params_.sublist("Printing");
    }
  
  TEST_FOR_EXCEPTION(linSolver_.ptr().get()==0, runtime_error,
                     "null linear solver object in NOXSolver ctor");

  TEST_FOR_EXCEPTION(statusTest_.get()==0, runtime_error,
                     "null status test object in NOXSolver ctor");

}




NOX::StatusTest::StatusType 
NOXSolver::solve(const NonlinearOperator<double>& F, 
                 TSFExtended::Vector<double>& solnVec) const 
{
  TimeMonitor timer(noxSolverTimer());

  Vector<double> x0 = F.getInitialGuess();
  RCP<NOX::TSF::Group> grp = rcp(new NOX::TSF::Group(x0, F, linSolver_));
  RCP<Teuchos::ParameterList> noxParams 
    = Teuchos::rcp(&params_, false);
  RCP<NOX::Solver::Generic> solver 
    = NOX::Solver::buildSolver(grp, statusTest_, noxParams);

  NOX::StatusTest::StatusType rtn = solver->solve();

  const NOX::TSF::Group* solnGrp 
    = dynamic_cast<const NOX::TSF::Group*>(&(solver->getSolutionGroup()));

  TEST_FOR_EXCEPTION(solnGrp==0, runtime_error,
                     "Solution group could not be cast to NOX::TSF::Group");

  const NOX::TSF::Vector* x 
    = dynamic_cast<const NOX::TSF::Vector*>(&(solnGrp->getX()));

  TEST_FOR_EXCEPTION(x==0, runtime_error,
    "Solution vector could not be cast to NOX::TSF::Vector");
  
  solnVec = x->getTSFVector();

  return rtn;
}
