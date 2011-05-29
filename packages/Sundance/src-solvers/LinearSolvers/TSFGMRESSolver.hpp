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

#ifndef TSFGMRESSOLVER_HPP
#define TSFGMRESSOLVER_HPP

#include "SundanceDefs.hpp"
#include "TSFKrylovSolver.hpp"
#include "SundanceHandleable.hpp"
#include "SundancePrintable.hpp"
#include "Teuchos_Describable.hpp"
#include "TSFLinearCombinationDecl.hpp"

namespace TSFExtended
{
  using namespace Teuchos;

  /**
   *
   */
  template <class Scalar>
  class GMRESSolver : public KrylovSolver<Scalar>,
                      public Sundance::Handleable<LinearSolverBase<Scalar> >,
                      public Printable,
                      public Describable
  {
  public:
    /** */
    GMRESSolver(const ParameterList& params = ParameterList())
      : KrylovSolver<Scalar>(params) {;}
    /** */
    GMRESSolver(const ParameterList& params,
                const PreconditionerFactory<Scalar>& precond)
      : KrylovSolver<Scalar>(params, precond) {;}

    /** */
    virtual ~GMRESSolver(){;}

    /** \name Printable interface */
    //@{
    /** Write to a stream  */
    void print(std::ostream& os) const 
    {
      os << description() << "[" << std::endl;
      os << this->parameters() << std::endl;
      os << "]" << std::endl;
    }
    //@}

    /** */
    int getKSpace() const {return this->parameters().template get<int>("Restart");}
    
    /** \name Describable interface */
    //@{
    /** Write a brief description */
    std::string description() const {return "GMRESSolver";}
    //@}

    /** \name Handleable interface */
    //@{
    /** Return a ref count pointer to a newly created object */
    virtual RCP<LinearSolverBase<Scalar> > getRcp() 
    {return rcp(this);}
    //@}
    
  protected:

    /** */
    virtual SolverState<Scalar> solveUnprec(const LinearOperator<Scalar>& A,
                                            const Vector<Scalar>& rhs,
                                            Vector<Scalar>& soln) const ;

    
  };

  template <class Scalar> inline
  SolverState<Scalar> GMRESSolver<Scalar>
  ::solveUnprec(const LinearOperator<Scalar>& A,
                const Vector<Scalar>& b,
                Vector<Scalar>& soln) const
  {
    int myRank = MPIComm::world().getRank();

    int maxiters = this->getMaxiters();
    int kSpace = getKSpace();
    Scalar tol = this->getTol();
    int verbosity = this->verb();

    if (verbosity > 1)
      {
        std::cerr << "GMRES solver" << std::endl;
        std::cerr << "Max iterations " << maxiters << std::endl;
        std::cerr << "Krylov subspace size " << kSpace<< std::endl;
        std::cerr << "Convergence tolerance " << tol << std::endl;
      }

    // following GMRES from Matlab
    Scalar normOfB = sqrt(b.dot(b));

    /* check for trivial case of zero rhs */
    if (normOfB < tol) 
      {
        soln = b.space().createMember();
        soln.zero();
        return SolverState<Scalar>(SolveConverged, 
                                   "yippee!!", 0, 0.0); 
      }

    soln = b.copy();


    Vector<Scalar> x0 = b.copy();
    Vector<Scalar> r0 = b.space().createMember(); 
    Vector<Scalar> vh = b.space().createMember(); 
    Vector<Scalar> tmp = b.space().createMember(); 
    Vector<Scalar> residVec = b.space().createMember(); 
    Vector<Scalar> u = b.space().createMember(); 
    Vector<Scalar> vrf = b.space().createMember(); 

    std::vector<Scalar> h(kSpace+1);
    std::vector<Scalar> f(kSpace+1);
    std::vector<Scalar> q(kSpace+1);
    std::vector<Scalar> mtmp(kSpace+1);
    std::vector<Scalar> y(kSpace+1);

    std::vector<Vector<Scalar> > V(kSpace + 1);
    std::vector<Vector<Scalar> > W(kSpace + 1);
    std::vector<std::vector<Scalar> > QT(kSpace + 1);
    std::vector<std::vector<Scalar> > R(kSpace + 1);



    for (int k = 0; k < V.size(); k++) 
      {		
        V[k] = A.domain().createMember(); // V = n x (m+1)
        W[k] = A.domain().createMember(); // W = n x (m+1)
        QT[k].resize(kSpace+1);
        R[k].resize(kSpace+1); // R = (m+1)x(m+1)
      }


    Scalar relTol = tol * normOfB; // relative tolerance
    Scalar phibar;
    Scalar rt;
    Scalar c;
    Scalar s;
    Scalar temp;

    int CONV = 0; // not converged yet

    // r0 =  b - A*x0;
    r0 = b - A*x0;
    residVec = r0.copy();
    Scalar normOfResidVec;
    normOfResidVec = residVec.norm2();
  
    // Outer loop i = 1 : maxIters unless convergence (or failure)
    //  for (int i=0; i<maxIters_; i++)
    int iter = 0;
    //    int i = 0;
    int j = 0;
    while (iter < maxiters)
      {
        for (int z=0; z<h.size(); z++) 
          {
            h[z]=0.0;
            f[z]=0.0;
          }
        for (int z = 0; z < V.size(); z++) 
          {
            V[z].zero();
            W[z].zero();
            for (int zz=0; zz<QT[z].size(); zz++) 
              {QT[z][zz]=0.0; R[z][zz]=0.0;}
          }

        vh = residVec.copy();
        h[0] = vh.norm2();
        double newtemp = 1.0 / h[0];
        V[0] = newtemp*vh;
        QT[0][0] = 1.0;
        phibar = h[0];

        // inner loop from 0:restart
        // for(int j=0; j<kSpace; j++)
        j = 0;
        while ((j < kSpace) & (iter < maxiters))
          {
            u = A*V[j];
            for( int k=0; k<=j; k++)
              {
                h[k] = V[k].dot(u);
                u = u - h[k] * V[k];
              }
					
            h[j+1] = u.norm2();
            V[j+1] = (1.0 / h[j+1]) * u;

            for (int k=0; k<=j; k++)
              {
                for (int z=0; z<=j; z++) 
                  R[k][j] = R[k][j] + QT[k][z] * h[z];
              }

            rt = R[j][j];

            // find cos(theta) and sin(theta) of Givens rotation
            if (h[j+1] == 0)
              {
                c = 1.0; // theta = 0 
                s = 0.0;
              }
            else if (fabs(h[j+1]) > fabs(rt))
              {
                temp = rt / h[j+1];
                // pi/4 < theta < 3pi/4
                s = 1.0 / sqrt(1.0 + fabs(temp)*fabs(temp)); 
                c = - temp * s;
              }
            else
              {
                temp = h[j+1] / rt;
                // -pi/4 <= theta < 0 < theta <= pi/4
                c = 1.0 / sqrt(1.0 + fabs(temp)*fabs(temp)); 
                s = - temp * c;
              }

            R[j][j] = c * rt - s * h[j+1]; // left out conj on c and s


            for(int k=0; k<=j; k++)
              q[k] = QT[j][k];
            for(int k=0; k<=j; k++)
              {
                QT[j][k] = c * q[k];
                QT[j+1][k] = s * q[k];
              }
            QT[j][j+1] = -s;
            QT[j+1][j+1] =  c;
            f[j] = c * phibar;
            phibar = s * phibar;

            if (j < kSpace-1)
              {
                W[j] = V[j].copy();
                for(int k=0; k<=j-1; k++)
                  W[j] = W[j] - R[k][j]*W[k];
                W[j] = (1.0 / R[j][j]) * W[j];

                soln = soln + f[j]*W[j];
              }
            else
              {
                for (int zz=0; zz<mtmp.size(); zz++) mtmp[zz]=0.0;
                // back solve to get tmp vector to form vrf
                mtmp[j] = f[j] / R[j][j];
                for(int k=j-1; k>=0; k--)
                  {
                    mtmp[k] = f[k];
                    for(int z=k+1; z<=j; z++)
                      mtmp[k] = mtmp[k] - R[k][z]*mtmp[z];
                    mtmp[k] = mtmp[k] / R[k][k];
                  }
							
                vrf.zero();
                for(int k=0; k<=j; k++)
                  vrf = vrf + mtmp[k] * V[k];
                soln = x0 + vrf;
              }
					
            // update current resid norm
            tmp.zero();
            tmp = A*soln;
            residVec = b - tmp;
            normOfResidVec = residVec.norm2();
            

            if (myRank==0 && verbosity > 1 ) 
              {
                std::cerr << "GMRES: iteration=";
                std::cerr.width(8);
                std::cerr << iter;
                std::cerr.width(20);
                std::cerr << "scaled resid=" << normOfResidVec/normOfB << std::endl;
              }

            // check for convergence
            if (normOfResidVec < relTol)
              {
                if (j < kSpace-1)
                  {
                    // compute more accurate soln to test convergence
                    for (int zz=0; zz<y.size(); zz++) y[zz]=0.0;

                    // back solve to get y(0:j) = R(0:j,0:j) \ f(0:j);
                    y[j] = f[j] / R[j][j];
                    for(int k=j-1; k>=0; k--)
                      {
                        y[k] = f[k];
                        for(int z=k+1; z<=j; z++)
                          y[k] = y[k] - R[k][z]*y[z];
                        y[k] = y[k] / R[k][k];
                      }
                    // soln = x0 + V(:,0:j) * y(0:j)
                    soln = x0.copy();
                    for(int k=0; k<=j; k++)
                      soln = soln + y[k] * V[k];
									
                    tmp.zero();
                    tmp = A*soln;
                    residVec = b - tmp;
                    normOfResidVec = residVec.norm2();
                  }
                // test for convergence
                if (normOfResidVec < relTol)
                  {
                    // we're done
                    CONV = 1; // converged
                    if (verbosity > 0 && myRank==0)
                      {
                        std::cerr << "GMRES converged (1) in " << iter+1
                             << " iters: final scaled resid = " 
                             <<  normOfResidVec/normOfB << std::endl;
                      }
                    SolverState<Scalar> rtn(SolveConverged, 
                                            "yippee!!", iter+1, 
                                            normOfResidVec/normOfB);
                    return rtn;
                  }
              }
            
            j++;
            iter++;
          } // end inner loop
			
        if (CONV)
          {
            if (verbosity > 0 && myRank==0)
              {
                std::cerr << "GMRES converged (1) in " << iter+1
                     << " iters: final scaled resid = " 
                     <<  normOfResidVec/normOfB << std::endl;
              }
            SolverState<Scalar> rtn(SolveConverged, 
                                    "yippee!!", iter+1, 
                                    normOfResidVec/normOfB);
            return rtn;
          }


        else if (!CONV & (iter < maxiters))
          {
            // not converged yet; update x0 and resid and restart
            x0 = soln.copy();
            tmp.zero();
            tmp = A*x0;
            residVec = b - tmp;
            normOfResidVec = residVec.norm2();
            
            if (verbosity > 1 && myRank==0)
              {
                std::cerr << "GMRES restarting: current scaled resid = "
                     << normOfResidVec/normOfB << std::endl;
              }
          }
			
      } // end outer loop
    
    SolverState<Scalar> rtn(SolveFailedToConverge, 
                            "GMRES failed to converge", 
                            maxiters, normOfResidVec/normOfB);
    return rtn;
	}
}

#endif
