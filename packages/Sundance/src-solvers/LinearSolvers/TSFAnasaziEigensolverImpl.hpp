/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// **********************************************************************/
/* @HEADER@ */

#ifndef TSFANASAZIEIGENSOLVER_IMPL_HPP
#define TSFANASAZIEIGENSOLVER_IMPL_HPP

#include "SundanceDefs.hpp"
#include "TSFAnasaziEigensolverDecl.hpp" 
#include "TSFParameterListPreconditionerFactory.hpp" 
#include "TSFPreconditionerFactory.hpp" 
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziSimpleLOBPCGSolMgr.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziThyraAdapter.hpp"
#include "TSFAnasaziAdapter.hpp"


namespace TSFExtended
{
using Teuchos::ParameterList;
using Anasazi::SimpleMV;
/** */
template <class MV, class OP> 
class InitTraits
{
public:
  /** */
  static RCP<OP> opPtr(const LinearOperator<double>& A);

  /** */
  static RCP<MV> makeMV(int numVecs, const VectorSpace<double>& space);

  /** */
  static Vector<double> vec(const RCP<MV>& mv, int i);
};


/** */
template <> class InitTraits<SimpleMV, LinearOperator<double> >
{
public:
  typedef SimpleMV            MV;
  typedef LinearOperator<double>            OP;

  /** */
  static RCP<OP> opPtr(const LinearOperator<double>& A)
    {
      if (A.ptr().get() != 0)
        return rcp(new LinearOperator<double>(A));
      else
      {
        RCP<LinearOperator<double> > rtn;
        return rtn;
      }
    }

  /** */
  static RCP<MV> makeMV(int blockSize, const VectorSpace<double>& space)
    {
      RCP<MV> mv = rcp(new MV(blockSize));
      for (int i=0; i<blockSize; i++) (*mv)[i] = space.createMember();
      return mv;
    }

  /** */
  static Vector<double> vec(const RCP<MV>& mv, int i)
    {
      return (*mv)[i];
    }

  
};


/** */
template <> class InitTraits<MultiVectorBase<double>, LinearOpBase<double> >
{
public:
  typedef Thyra::MultiVectorBase<double>         MV;
  typedef Thyra::LinearOpBase<double>            OP;

  /** */
  static RCP<OP> opPtr(const LinearOperator<double>& A)
    {
      return A.ptr();
    }

  /** */
  static RCP<MV> makeMV(int blockSize, const VectorSpace<double>& space)
    {
      RCP<const Thyra::VectorSpaceBase<double> > mvSpace = space.ptr();
      return Thyra::createMembers( *mvSpace, blockSize );
    }

  /** */
  static Vector<double> vec(const RCP<MV>& mv, int i)
    {
      return mv->col(i);
    }
};





template <class Scalar>  
inline void AnasaziEigensolver<Scalar>::solve(
  const LinearOperator<Scalar>& K,
  const LinearOperator<Scalar>& M,
  Array<Vector<Scalar> >& evecs,
  Array<std::complex<Scalar> >& ew) const 
{
//#define USE_THYRA_MV

#ifdef USE_THYRA_MV
  typedef Thyra::MultiVectorBase<Scalar>         MV;
  typedef Thyra::LinearOpBase<Scalar>            OP;
#else
  typedef SimpleMV            MV;
  typedef LinearOperator<Scalar>            OP;
#endif
  typedef Anasazi::MultiVecTraits<Scalar,MV>     MVT;
  typedef Anasazi::OperatorTraits<Scalar,MV,OP>  OPT;

  TimeMonitor timer(solveTimer());
  VectorSpace<Scalar> KDomain = K.domain();

  RCP<OP> KPtr = InitTraits<MV, OP>::opPtr(K);
  RCP<OP> MPtr = InitTraits<MV, OP>::opPtr(M);



  
  // Eigensolver parameters
  std::string method = this->params().get<string>("Method");
  int numEigs = this->params().get<int>("Number of Eigenvalues");
  int blockSize = this->params().get<int>("Block Size");
  bool usePrec = this->params().get<bool>("Use Preconditioner");
  bool hermitian = this->params().get<bool>("Is Hermitian");


  
  /* Make a multivector with row space = domain of K, column 
   * space = multiVec Space*/
  RCP<MV> mv = InitTraits<MV, OP>::makeMV(blockSize, KDomain);

  /* Fill the multivector with random values */
  MVT::MvRandom( *mv );

  /* Create a preconditioner */
  ParameterList precParams = this->params().sublist("Preconditioner");
  PreconditionerFactory<double> precFactory 
    = new ParameterListPreconditionerFactory(precParams);

  LinearOperator<Scalar> P;
  if (usePrec) 
  {
    TimeMonitor pTimer(precondBuildTimer());
    P = precFactory.createPreconditioner(K).right();
  }

  /* Create eigenproblem */
  RCP<Anasazi::Eigenproblem<Scalar,MV,OP> > problem;

  if (MPtr.get() != 0)
  {
    problem =
      rcp( new Anasazi::BasicEigenproblem<Scalar,MV,OP>(KPtr, MPtr, mv) );
  }
  else
  {
    problem =
      rcp( new Anasazi::BasicEigenproblem<Scalar,MV,OP>(KPtr, mv) );
  }

  ParameterList eigParams = this->params();
  problem->setHermitian(hermitian);
  problem->setNEV(numEigs);
  if (usePrec) problem->setPrec(InitTraits<MV, OP>::opPtr(P));

  bool ret = problem->setProblem();
  TEST_FOR_EXCEPTION(!ret, std::runtime_error,
    "Eigenproblem not setup correctly");
  

  // Create the solver manager
  RCP<Anasazi::SolverManager<Scalar,MV,OP> > MySolverMan;
  if (method=="Block Davidson")
  {
    MySolverMan = rcp(new Anasazi::BlockDavidsonSolMgr<Scalar,MV,OP>(problem, eigParams));
  }
  else if (method=="Block Krylov Schur")
  {
    MySolverMan = rcp(new Anasazi::BlockKrylovSchurSolMgr<Scalar,MV,OP>(problem, eigParams));
  }
  else if (method=="Simple LOBPCG")
  {
    MySolverMan = rcp(new Anasazi::SimpleLOBPCGSolMgr<Scalar,MV,OP>(problem, eigParams));
  }
  else if (method=="LOBPCG")
  {
    MySolverMan = rcp(new Anasazi::LOBPCGSolMgr<Scalar,MV,OP>(problem, eigParams));
  }
  else
  {
    TEST_FOR_EXCEPTION(true, std::runtime_error,
      "solver method [" << method << "] not recognized");
  }

  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMan->solve();
  Out::os() << "return code = " << returnCode << endl;
  TEST_FOR_EXCEPTION(returnCode != Anasazi::Converged, 
    std::runtime_error, "Anasazi did not converge!");
  
  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<Scalar,MV> sol = problem->getSolution();
  RCP<MV> evecs_mv = sol.Evecs;
  int numev = sol.numVecs;
  
  /* Copy the columns of the eigenvector MV into an array of TSF vectors */
  ew.resize(numev);
  evecs.resize(numev);

  for (int i=0; i<numev; i++)
  {
    Vector<Scalar> tmp = InitTraits<MV, OP>::vec(evecs_mv, i);

    evecs[i] = KDomain.createMember();
    evecs[i].acceptCopyOf(tmp);
    /* record the associated eigenvalue. The matrix is Hermitian so
     * we know the eigenvalue is real. */
    //evals[i] = sol.Evals[i].realpart;
    // if matrix might not be hermitian
    ew[i].real() = sol.Evals[i].realpart;
    ew[i].imag() = sol.Evals[i].imagpart;
  }
}


}


#endif
