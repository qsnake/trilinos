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

#include "TSFPoissonBoltzmannOp.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorImpl.hpp"
#include "TSFLinearOperatorImpl.hpp"
#endif

using namespace TSFExtended;
using namespace Teuchos;


PoissonBoltzmannOp::PoissonBoltzmannOp(int nLocal, const VectorType<double>& vecType)
  : NonlinearOperatorBase<double>(), J_(nLocal, vecType), importer_(),
    uLeftBC_(0.0), uRightBC_(2.0*log(cosh(1.0/sqrt(2.0))))
{
  setDomainAndRange(J_.domain(), J_.range());

  int rank = MPIComm::world().getRank();
  int nProc = MPIComm::world().getNProc();
  if (nProc > 1)
    {
      Array<int> ghosts;
      int low = J_.domain().lowestLocallyOwnedIndex();
      int high = low + J_.domain().numLocalElements();
      if (rank != nProc - 1)
        {
          ghosts.append(high);
        }
      if (rank != 0) 
        {
          ghosts.append(low-1);
        }

      importer_ = vecType.createGhostImporter(J_.domain(), ghosts.size(), &(ghosts[0]));
    }
  else
    {
      importer_ = vecType.createGhostImporter(J_.domain(), 0, 0);
    }
}

Vector<double> PoissonBoltzmannOp::getInitialGuess() const
{
  Vector<double> rtn = J_.domain().createMember();

  rtn.setToConstant(0.5);

  return rtn;
}


LinearOperator<double> 
PoissonBoltzmannOp::computeJacobianAndFunction(Vector<double>& functionValue) const 
{
  J_.setEvalPoint(currentEvalPt());

  RCP<GhostView<double> > u;
  importer_->importView(currentEvalPt(), u);

  int low = J_.domain().lowestLocallyOwnedIndex();
  int high = low + J_.domain().numLocalElements();

  functionValue = J_.range().createMember();
  double h= J_.h();

  for (int r=low; r<high; r++)
    {
      double u_i = u->getElement(r);
      double f = 0.0;
      if (r==0) 
        {
          f = u_i - uLeftBC_;
        }
      else if (r==J_.domain().dim()-1)
        {
          f = u_i - uRightBC_;
        }
      else
        {
          double u_plus = u->getElement(r+1);
          double u_minus = u->getElement(r-1);
          f = (u_plus + u_minus - 2.0*u_i)/h/h - exp(-u_i);
        }
      functionValue.setElement(r, f);
    }

  return J_.getOp();
}

Vector<double> PoissonBoltzmannOp::exactSoln() const
{
  Vector<double> rtn = J_.domain().createMember();

  int low = J_.domain().lowestLocallyOwnedIndex();
  int high = low + J_.domain().numLocalElements();
  double h= J_.h();
  
  for (int r=low; r<high; r++)
    {
      double x = r*h;
      double u = 2.0*log(cosh(x/sqrt(2.0)));
      std::cerr << x << " " << u << std::endl;
      rtn.setElement(r, u);
    }

  return rtn;
}

