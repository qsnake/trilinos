//@HEADER
// ***********************************************************************
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
// ***********************************************************************
//@HEADER

#include "Teuchos_GlobalMPISession.hpp"
#include "TSFVectorDecl.hpp"
#include "TSFLinearCombinationDecl.hpp"
#include "TSFVectorType.hpp"
#include "TSFVectorSpaceDecl.hpp"
#include "TSFEpetraVectorType.hpp"
#include "Teuchos_Time.hpp"
#include "TSFHeatOperator1D.hpp"
#include "TSF_NVector.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#include "TSFLinearSolverImpl.hpp"
#include "TSFLinearCombinationImpl.hpp"
#include "TSFCommonOperatorsDecl.hpp"
#endif



using namespace Teuchos;
using namespace TSFExtended;
using namespace TSFExtendedOps;


#ifdef HAVE_SUNDIALS
#include "ida.h"
#include "idaspgmr.h"
//#include "ida_impl.h"
#endif

class IDAHeat
{
public:
  /** */
  IDAHeat(int nLocal, const VectorType<double>& vecType 
          /*,
            const LinearSolver<double>& solver*/) ;

  /** */
  VectorSpace<double> space() const {return vecSpace_;}

  /** */
  void computeResidual(const double& time, 
                       const Vector<double>& u,
                       const Vector<double>& uDot,
                       Vector<double>& resid) const ;

  /** */
  LinearOperator<double> jacobian(const double& time, 
                                  const double& c_j,
                                  const Vector<double>& u,
                                  const Vector<double>& uDot) const ;

  /** */
  void initialConditions(Vector<double>& u,
                         Vector<double>& uDot,
                         Vector<double>& resid) const ;

  /** */
  SolverState<double> solve(const double& t,
                            const double& cj, 
                            const Vector<double>& u,
                            const Vector<double>& uDot,
                            const Vector<double>& b,
                            Vector<double>& x) const ;

  /** */
  Vector<double> exactSoln(const double& t) const ;

private:
  RCP<GhostImporter<double> > importer_;
  double h_;
  VectorType<double> vecType_;
  VectorSpace<double> vecSpace_;
  mutable HeatOperator1D opBuilder_;
  /*  LinearSolver<double> solver_; */
};

#ifdef HAVE_SUNDIALS

static int check_flag(void *flagvalue, char *funcname, int opt);

int resTrilinos(realtype tres, N_Vector uu, N_Vector up,
                N_Vector resval, void *rdata)
{
  IDAHeat* h = reinterpret_cast<IDAHeat*>(rdata);
  Vector<realtype> u = toTrilinos(uu);
  Vector<realtype> uDot = toTrilinos(up);
  Vector<realtype> resid = toTrilinos(resval);
  h->computeResidual(tres, u, uDot, resid);
  toTrilinos(resval) = resid;
  return 0;
}

#ifdef BLARF
int solveTrilinos(IDAMem idaMem, N_Vector b, N_Vector w,
                  N_Vector u, N_Vector uDot, N_Vector resid)
{
  double cj = idaMem->ida_cj;
  double tn = idaMem->ida_tn;
  void* rdata = idaMem->ida_rdata;
  IDAHeat* h = reinterpret_cast<IDAHeat*>(rdata);

  Vector<double> x; 
  SolverState<double> state = h->solve(t, cj, toTrilinos(u), toTrilinos(uDot),
                                       toTrilinos(b), x);
  toTrilinos(b).acceptCopyOf(x);
  return state.finalState() != SolveConverged;
}
#endif

#endif

IDAHeat::IDAHeat(int nLocal, const VectorType<double>& vecType 
        /*,
          const LinearSolver<double>& solver*/) 
  : importer_(),
    h_(0.0), 
    vecType_(vecType),
    vecSpace_(),
    opBuilder_(nLocal, vecType) /*,
                                  solver_(solver)*/
{
  const double pi = 4.0*atan(1.0);
  vecSpace_ = opBuilder_.domain();
  int n = vecSpace_.dim();
  h_ = pi/(n - 1.0);
  
  int rank = MPIComm::world().getRank();
  int nProc = MPIComm::world().getNProc();
  int low = vecSpace_.lowestLocallyOwnedIndex();

  Array<int> ghosts;  
  if (nProc > 1)
    {
      if (rank==0)
        {
          ghosts=tuple(nLocal);
        }
      else if (rank==nProc-1)
        {
          ghosts=tuple(low-1);
        }
      else
        {
          ghosts=tuple(low-1, low+nLocal);
        }
    }
  if (ghosts.size()>0)
    {
      importer_ = vecType_.createGhostImporter(vecSpace_, 
                                               ghosts.size(), 
                                               &(ghosts[0]));
    }
  else
    {
      importer_ = vecType_.createGhostImporter(vecSpace_, 
                                               ghosts.size(), 0);
    }
}

void IDAHeat::computeResidual(const double& time, 
                              const Vector<double>& u,
                              const Vector<double>& uDot,
                              Vector<double>& resid) const
{
  int low = vecSpace_.lowestLocallyOwnedIndex();
  int nLocal = vecSpace_.numLocalElements();
  int n = vecSpace_.dim();

  /* Initialize rr to uu, to take care of boundary equations. */
  resid.acceptCopyOf(u);

  RCP<GhostView<double> > gu;
  RCP<GhostView<double> > guDot;
  importer_->importView(u, gu);
  importer_->importView(uDot, guDot);

  for (int i=0; i<nLocal; i++)
    {
      int g = i + low;
      if (g==0 || g==(n-1)) 
        {
          resid.setElement(g, u.getElement(g));
        }
      else
        {
          double dif1 = gu->getElement(g-1) + gu->getElement(g+1)
            - 2.0*gu->getElement(g);
          double MuDot = (guDot->getElement(g-1) + guDot->getElement(g+1)
                          + 4.0*guDot->getElement(g))/6.0;
          double r = MuDot - dif1/h_/h_;
          resid.setElement(g, r);
        }
    }
}


void IDAHeat::initialConditions(Vector<double>& u,
                                Vector<double>& uDot,
                                Vector<double>& resid) const
{
  
  int rank = MPIComm::world().getRank();
  int nProc = MPIComm::world().getNProc();

  int low = vecSpace_.lowestLocallyOwnedIndex();
  int nLocal = vecSpace_.numLocalElements();
  int n = vecSpace_.dim();

  /* Initialize u on all grid points. */
  for (int i=0; i<nLocal; i++) 
    {
      int g = i + low;
      double uVal = sin(g*h_);
      u.setElement(g, uVal);
    }

  /* resHeat sets res to negative of ODE RHS values at interior points. */
  uDot.setToConstant(0.0);
  computeResidual(0.0, u, uDot, resid);
  uDot.acceptCopyOf(-1.0*resid);


  /* Set uDot to zero at boundary points. */
  if (rank==0) uDot.setElement(0, 0.0);
  if (rank==nProc-1) uDot.setElement(n-1, 0.0);
}

LinearOperator<double> IDAHeat::jacobian(const double& /* time */, 
                                         const double& c_j,
                                         const Vector<double>& /* u */,
                                         const Vector<double>& /* uDot */) const
{
  opBuilder_.set_cj(c_j);
  return opBuilder_.getOp();
}

/*
SolverState<double> IDAHeat::solve(const double& t,
                                   const double& cj, 
                                   const Vector<double>& u,
                                   const Vector<double>& uDot,
                                   const Vector<double>& b,
                                   Vector<double>& x) const
{
  LinearOperator<double> J = jacobian(t, cj, u, uDot); 
  return solver_.solve(J, b, x);
}
*/


Vector<double> IDAHeat::exactSoln(const double& t) const
{
  Vector<double> rtn = vecSpace_.createMember();

  int low = vecSpace_.lowestLocallyOwnedIndex();
  int nLocal = vecSpace_.numLocalElements();


  for (int i=0; i<nLocal; i++) 
    {
      int g = i + low;
      double uVal = sin(g*h_)*exp(-t);
      rtn.setElement(g, uVal);
    }
  return rtn;
}

int main(int argc, char *argv[]) 
{
  try
    {
      GlobalMPISession session(&argc, &argv);

      int n = 201;
      


      int rank = MPIComm::world().getRank();
      int nProc = MPIComm::world().getNProc();
#ifndef HAVE_SUNDIALS
      cout << "sundials not present... test INACTIVE" << std::endl;
#else

      IDAHeat model(n, new EpetraVectorType());
      VectorSpace<double> space = model.space();

      N_Vector uu = N_VNew_Trilinos(space);
      N_Vector up = N_VNew_Trilinos(space);
      N_Vector res = N_VNew_Trilinos(space);
      N_Vector constraints = N_VNew_Trilinos(space);

      toTrilinos(constraints).setToConstant(1.0);

      model.initialConditions(toTrilinos(uu),
                              toTrilinos(up),
                              toTrilinos(res));

      double t0 = 0.0;
      double t1 = 0.1;
      double rtol = 1.0e-6;
      double atol = 1.0e-6;



      /* Call IDACreate and IDAMalloc to initialize solution */

      void* mem = IDACreate();
      TEST_FOR_EXCEPTION(check_flag((void*) mem, "IDACreate", 0), 
                         std::runtime_error,
                         "IDACreate returned null object");

      int ier = IDASetRdata(mem, &model);
      TEST_FOR_EXCEPTION(check_flag(&ier, "IDASetRdata", 1), std::runtime_error,
                         "IDASetRdata returned error code " << ier);

  //     ier = IDASetConstraints(mem, constraints);
//       TEST_FOR_EXCEPTION(check_flag(&ier, "IDASetConstraints", 1), 
//                          std::runtime_error,
//                          "IDASetConstraints returned error code " << ier);

      N_VDestroy_Trilinos(constraints);

      ier = IDAMalloc(mem, resTrilinos, t0, uu, up, IDA_SS, rtol, &atol);
      TEST_FOR_EXCEPTION(check_flag(&ier, "IDAMalloc", 1), std::runtime_error,
                         "IDAMalloc returned error code " << ier);

      /* Call IDASpgmr to specify the linear solver. */
      /* The solver is unpreconditioned GMRES */

      int krylov = 1000;
      ier = IDASpgmr(mem, krylov);
      TEST_FOR_EXCEPTION(check_flag(&ier, "IDASpgmr", 1), std::runtime_error,
                         "IDASpgmr returned error code " << ier);

      

      /* Loop over output times, call IDASolve, and print results. */
      
      int nOut = 11;
      double tout=0.0;
      double uMax = toTrilinos(uu).normInf();
      double uErr = (model.exactSoln(tout)-toTrilinos(uu)).normInf();
      double maxErr = uErr;
      if (rank==0) printf("%f  %14.7g  %14.7g\n", tout, 
                          uMax, uErr);
      double tret;
      int iout;
      for (tout = t1, iout = 1; 
           iout <= nOut ; 
           iout++, tout += 0.1) 
        {
          ier = IDASolve(mem, tout, &tret, uu, up, IDA_NORMAL);
          TEST_FOR_EXCEPTION(check_flag(&ier, "IDASolve", 1), std::runtime_error,
                             "IDASolve returned error code " << ier);

          uMax = toTrilinos(uu).normInf();
          uErr = (model.exactSoln(tout)-toTrilinos(uu)).normInf();
          if (rank==0) printf("%f  %14.7g  %14.7g\n", tout, 
                              uMax, uErr);
          if (uErr > maxErr) maxErr = uErr;
      }

      /* Free Memory */
      
      IDAFree(mem);
      
      N_VDestroy_Trilinos(uu);
      N_VDestroy_Trilinos(up);
      N_VDestroy_Trilinos(res);

      double tol = 1.0e-5;
      if (maxErr > tol)
        {
          std::cerr << "IDA heat eqn test FAILED" << std::endl;
        }
      else
        {
          std::cerr << "IDA heat eqn test PASSED" << std::endl;
        }

#endif 
    }
  catch(std::exception& e)
    {
      cout << "Caught exception: " << e.what() << std::endl;
    }
}



#ifdef HAVE_SUNDIALS
/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr,
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1);
  } else if (opt == 1) {
    /* Check if flag < 0 */
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr,
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return(1);
    }
  } else if (opt == 2 && flagvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1);
  }

  return(0);
}
#endif
