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

#ifndef SUNDANCE_LINEAREIGENPROBLEM_H
#define SUNDANCE_LINEAREIGENPROBLEM_H

#include "SundanceDefs.hpp"
#include "SundanceLinearProblem.hpp"
#include "SundanceEigensolution.hpp"
#include "SundanceFunctionalEvaluator.hpp"
#include "TSFEigensolver.hpp"

namespace Sundance
{
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using namespace TSFExtended;


/**
 *
 */
class LinearEigenproblem
{
public:
  /** */
  LinearEigenproblem(){}
  /** */
  LinearEigenproblem(
    const Mesh& mesh, const Expr& eqn,
    const Expr& v, const Expr& u,
    const VectorType<double>& vecType) ;
  /** */
  LinearEigenproblem(
    const Mesh& mesh, const Expr& eqn,
    const Expr& v, const Expr& u,
    const VectorType<double>& vecType,
    bool lumpMass) ;
  /** */
  LinearEigenproblem(
    const Mesh& mesh, const Expr& eqn,
    const Expr& massExpr,
    const Expr& v, const Expr& u,
    const VectorType<double>& vecType,
    bool lumpMass) ;

  /** */
  Eigensolution solve(const Eigensolver<double>& solver) const ;

  /** */
  LinearOperator<double> getK() const {return kProb_.getOperator();}

  /** */
  LinearOperator<double> getM() const {return mProb_.getOperator();}
    
  
private:
  /** */
  Array<Expr> makeEigenfunctions(Array<Vector<double> >& ev) const ;

  /** */
  LinearProblem makeMassProb(
    const Mesh& mesh,
    const Expr& massExpr,
    const Expr& v, const Expr& u,
    const VectorType<double>& vecType) const ;

  /** */
  LinearOperator<double> 
  lumpedOperator(const LinearOperator<double>& M) const ;

  bool lumpMass_;
  LinearProblem kProb_;
  LinearProblem mProb_;
  LinearOperator<double> M_;
  LinearOperator<double> MUnlumped_;
  DiscreteSpace discSpace_;
};




}


#endif
