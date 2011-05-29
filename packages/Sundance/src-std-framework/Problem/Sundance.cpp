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
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_MPIComm.hpp"
#include "Teuchos_MPISession.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceQuadratureIntegral.hpp"
#include "SundanceQuadratureEvalMediator.hpp"
#include "SundanceDiscreteFuncEvaluator.hpp"
#include "SundanceRefIntegral.hpp"
#include "SundancePathUtils.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceCellJacobianBatch.hpp"
#include "SundanceUnaryFunctor.hpp"
#include "SundanceGrouperBase.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceDefaultPath.hpp"
#include "SundanceDOFMapBuilder.hpp"
#include "SundanceVersionString.hpp"
#include "SundanceProductTransformation.hpp"
#include <unistd.h>
#ifndef _MSC_VER
#include <sys/unistd.h>
#else
#include "winprocess.h"
#endif
#ifdef HAVE_MPI
#include "mpi.h"
#endif

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#endif

namespace Sundance
{

void handleException(std::exception& e)
{SundanceGlobal::handleException(e);}


bool passFailTest(bool pass)
{return SundanceGlobal::passFailTest(pass);}

bool passFailTest(double error, double tol)
{return SundanceGlobal::passFailTest(error, tol);}


bool passFailTest(const std::string& statusMsg,
  bool status, double error, double tol)
{return SundanceGlobal::passFailTest(statusMsg, status, error, tol);}

int& testStatus() 
{return SundanceGlobal::testStatus();}


CommandLineProcessor& clp() 
{return SundanceGlobal::clp();}


int init(int* argc, char*** argv)
{return SundanceGlobal::init(argc, argv);}


int finalize() {return SundanceGlobal::finalize();}



void setOption(const std::string& optionName,
  int& value,
  const std::string& helpMsg)
{
  SundanceGlobal::setOption(optionName, value, helpMsg);
}

void setOption(const std::string& optionName,
  double& value,
  const std::string& helpMsg)
{
  SundanceGlobal::setOption(optionName, value, helpMsg);
}

void setOption(const std::string& optionName,
  std::string& value,
  const std::string& helpMsg)
{
  SundanceGlobal::setOption(optionName, value, helpMsg);
}

void setOption(const std::string& optionTrueName,
  const std::string& optionFalseName,
  bool& value,
  const std::string& helpMsg)
{
  SundanceGlobal::setOption(optionTrueName,
    optionFalseName,
    value,
    helpMsg);
}
} // namespace Sundance

static Time& totalTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("total Sundance time"); 
  return *rtn;
}

int SundanceGlobal::init(int* argc, char*** argv)
{

  try
  {
    /* start up MPI. In a serial run, this will be a no-op */
    //      MPISession::init(argc, argv);
    globalMPISession(argc, (char***) argv);

    /* Start a stopwatch. It will be stopped upon a call to finalize() */
    totalTimer().start();

    Tabs tab;

    /* read standard command line flags */
    std::string configFilename = "";

    bool defaultFpCheck = false;
    bool debugWait = false;
    bool showVersion = false;
    bool showBanner = true;
    bool showTimings = false;
    bool cmdFpCheck = defaultFpCheck;
    int defaultWorkSetSize = 400;
    int cmdWorkSetSize = defaultWorkSetSize;

    Assembler::workSetSize() = defaultWorkSetSize;

    clp().setOption("config", &configFilename, "Configuration file");
    clp().setOption("fpcheck", "nofpcheck", &cmdFpCheck, 
      "Check results of math lib calls in expr evals");
    clp().setOption("version", "noversion", &showVersion, 
      "Show Sundance version number and exit");
    clp().setOption("banner", "nobanner", &showBanner, 
      "Show Sundance banner on root processor at start of run");
    clp().setOption("timings", "notimings", &showTimings, 
      "Show timings at end of run");

    clp().setOption("workset", &cmdWorkSetSize, 
      "Work set size");

      
    clp().setOption("debug", "nodebug", &debugWait, "Whether to attach a debugger to this process, holding until 'wait' is set to 0");


    clp().throwExceptions(false);
    clp().recogniseAllOptions(false);

    CommandLineProcessor::EParseCommandLineReturn rtn 
      = clp().parse(*argc, (char**) *argv);

    TEST_FOR_EXCEPTION(rtn != CommandLineProcessor::PARSE_SUCCESSFUL,
      RuntimeError,
      "Command-line parsing failed");


    if (showVersion)
    {
      if (MPIComm::world().getRank()==0)
      {
        cout << "Simulation built using Sundance version " 
             << VersionString::number() 
             << " (" << VersionString::date() << ")" << std::endl;
     
        cout << "Sundance is copyright (C) 2005 Sandia National Laboratories and is"
             << std::endl;
        cout << "licensed under the GNU Lesser General Public License, version 2.1" << std::endl;
        cout << tab << std::endl;
        exit(0);
      }
    }
    if (showBanner && MPIComm::world().getRank()==0)
    {
      cout << "Simulation built using Sundance version " 
           << VersionString::number() 
           << " (" << VersionString::date() << ")" << std::endl;
      
      cout << "Sundance is copyright" 
           << std::endl << " (C) 2005-2008 Sandia National Laboratories " 
           << std::endl
           << " (C) 2007-2008 Texas Tech University"
           << std::endl;
      cout << "and is licensed under the GNU Lesser General Public License, version 2.1" << std::endl;
      cout << tab << std::endl;
    }

    if (!showTimings) skipTimingOutput() = true;

    //      debugWait = true;
    if (debugWait)
    {
      int wait=1;
      int pid = getpid();
      std::string myCommandName=((char**)(*argv))[0];
      std::string debugCmd = "ddd --gdb -x ~/.gdbinit " + myCommandName 
        + " " + Teuchos::toString(pid) + " &";
      cout << "launching " << debugCmd << std::endl;
      system(debugCmd.c_str());
      while (wait) {;}
    }


    bool worksetSetOnCmdLine = cmdWorkSetSize != defaultWorkSetSize;
    if (worksetSetOnCmdLine)
    {
      Assembler::workSetSize() = (int) cmdWorkSetSize;
    }
  }
  catch(std::exception& e)
  {
    handleException(e);
    return 1;
  }
  return 0;
} 


bool& SundanceGlobal::showStartupMessage()
{
#ifdef TRILINOS_6
  static bool rtn=false; 
  return rtn;
#else
  return MPISession::showStartupMessage();
#endif
}



void SundanceGlobal::handleException(std::exception& e)
{
  cout << "Sundance detected exception: " << std::endl;
  cout << e.what() << std::endl;
  cout << "test FAILED" << std::endl;
  testStatus() = -1;
}


int SundanceGlobal::finalize()
{
  totalTimer().stop();

  try
  {
    Tabs tab;
    if (false && MPIComm::world().getRank()==0)
    {
      cout << tab << "eval vector flops: " << EvalVector::totalFlops() << std::endl;
      cout << tab << "quadrature flops: " << QuadratureIntegral::totalFlops() << std::endl;
      cout << tab << "ref integration flops: " 
           << RefIntegral::totalFlops() << std::endl;
      cout << tab << "cell jacobian batch flops: " << CellJacobianBatch::totalFlops() << std::endl;
      cout << tab << "quadrature eval mediator: " << QuadratureEvalMediator::totalFlops() << std::endl;
    }
    /* we may need to skip timing summaries because of a Trilinos 6.0.x bug */
    if (!skipTimingOutput()) TimeMonitor::summarize();
    //  MPISession::finalize();
  }
  catch(std::exception& e)
  {
    handleException(e);
    return 1;
  }
  return 0;
}






bool SundanceGlobal::checkTest(double error, double tol)
{
  int myFail = error > tol;
  int anyFail = myFail;
  MPIComm::world().allReduce((void*) &myFail, (void*) &anyFail, 1, MPIComm::INT,
    MPIComm::SUM);
  return (anyFail == 0);
}

bool SundanceGlobal:: passFailTest(double error, double tol)
{
  bool pass;
  if (MPIComm::world().getRank()==0)
  {
    cout << "error norm = " << error << std::endl;
    cout << "tolerance = " << tol << std::endl;
  }
  pass = checkTest(error, tol);
  if (MPIComm::world().getRank()==0)
  {
    if (pass)
    {
      cout << "test PASSED" << std::endl;
    }
    else
    {
      cout << "test FAILED" << std::endl;
    }
  }
  testStatus() = pass!=true;
  return pass;
}


bool SundanceGlobal:: passFailTest(const std::string& statusMsg,
  bool status, double error, double tol)
{
  bool pass;
  if (MPIComm::world().getRank()==0)
  {

    cout << statusMsg << ": ";
    if (status) cout << "true" << std::endl;
    else cout << "false" << std::endl;
    cout << "error norm = " << error << std::endl;
    cout << "tolerance = " << tol << std::endl;
  }
  pass = checkTest(error, tol);
  if (MPIComm::world().getRank()==0)
  {
    if (status && pass)
    {
      cout << "test PASSED" << std::endl;
    }
    else
    {
      cout << "test FAILED" << std::endl;
    }
  }
  testStatus() = pass!=true;
  return pass;
}

bool SundanceGlobal:: passFailTest(bool pass)
{
  if (MPIComm::world().getRank()==0)
  {
    if (pass)
    {
      cout << "test PASSED" << std::endl;
    }
    else
    {
      cout << "test FAILED" << std::endl;
    }
  }
  testStatus() = pass!=true;
  return pass;
}


void SundanceGlobal::setOption(const std::string& optionName,
  int& value,
  const std::string& helpMsg)
{
  clp().setOption(optionName.c_str(), &value, helpMsg.c_str());
}

void SundanceGlobal::setOption(const std::string& optionName,
  double& value,
  const std::string& helpMsg)
{
  clp().setOption(optionName.c_str(), &value, helpMsg.c_str());
}

void SundanceGlobal::setOption(const std::string& optionName,
  std::string& value,
  const std::string& helpMsg)
{
  clp().setOption(optionName.c_str(), &value, helpMsg.c_str());
}

void SundanceGlobal::setOption(const std::string& optionTrueName,
  const std::string& optionFalseName,
  bool& value,
  const std::string& helpMsg)
{
  clp().setOption(optionTrueName.c_str(),
    optionFalseName.c_str(),
    &value,
    helpMsg.c_str());
}


