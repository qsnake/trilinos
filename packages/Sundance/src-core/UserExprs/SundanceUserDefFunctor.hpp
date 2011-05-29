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

#ifndef SUNDANCE_USERDEFFUNCTOR_H
#define SUNDANCE_USERDEFFUNCTOR_H

#include "SundanceDefs.hpp"
#include "SundanceMap.hpp"
#include "SundanceMultiSet.hpp"
#include "SundanceExceptions.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"



namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;




class EvalVector;
class EvalManager;

/**
 * UserDefFunctor defines an interface for callbacks used to implement
 * user-defined nonlinear operators in the Sundance Expr system.
 */
class UserDefFunctor
{
public:
  /** ctor */
  UserDefFunctor(const std::string& name, int domainDim, int rangeDim) ;

  /** */
  virtual ~UserDefFunctor(){;}

  /** */
  const std::string& name(int elemIndex) const {return elemNames_[elemIndex];}

  /** */
  const std::string& name() const {return name_;}


  /** */
  virtual void evaluationCallback(int nPoints, int maxDiffOrder,
    const double** in,
    double** out) const = 0 ;

  /** */
  virtual void eval0(const Array<double>& in, double* outVals) const ;

  /**
   * Evaluate the expression and its derivative. The values should be put into
   * the outVals array. The derivatives should be put into the outDerivs array,
   * ordered with the domain index running fastest. That is, 
   * \f[
   * outDerivs[i*N_R + j] = \frac{\partial F_i}{\partial q_j}
   * \f]
   */
  virtual void eval1(const Array<double>& in, double* outVals, 
    double* outDerivs) const ;

    

  /** */
  int domainDim() const {return domainDim_;}

  /** */
  int rangeDim() const {return rangeDim_;}

  /** */
  virtual int maxOrder() const = 0 ;

  /** */
  void reset() const ;

private:
  const std::string name_;
  Array<string> elemNames_;
  const int domainDim_;
  const int rangeDim_;
};


}


#endif
