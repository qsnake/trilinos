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

#ifndef SUNDANCE_FUNCTIONALEVALUATOR_H
#define SUNDANCE_FUNCTIONALEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "TSFVectorImpl.hpp"
#include "TSFVectorType.hpp"

namespace Sundance
{

using namespace TSFExtended;
using namespace Teuchos;

class Mesh;
class Assembler;

/** 
 * 
 */
class FunctionalEvaluator 
{
public:
  /** */
  FunctionalEvaluator();

  /** */
  FunctionalEvaluator(const Mesh& mesh, 
    const Expr& integral);

  /** */
  FunctionalEvaluator(const Mesh& mesh, 
    const Expr& integral,
    const Expr& bcs,
    const Expr& var,
    const Expr& varEvalPts,
    const VectorType<double>& vectorType);

  /** */
  FunctionalEvaluator(const Mesh& mesh, 
    const Expr& integral,
    const Expr& bcs,
    const Expr& vars,
    const Expr& varEvalPts,
    const Expr& fields,
    const Expr& fieldValues,
    const VectorType<double>& vectorType);


  /** */
  double evaluate() const ;

  /** */
  Expr evalGradient(double& value) const ;

  /** */
  double fdGradientCheck(double h) const ;
          

private:

  /** */
  Vector<double> evalGradientVector(double& value) const ;
      
  /** */
  RCP<Assembler> assembler_;
      
  /** */
  mutable Expr varValues_;

  /** */
  VectorType<double> vecType_;
      
  /** */
  mutable Array<Vector<double> > gradient_;
      
};

/** */
double evaluateIntegral(const Mesh& mesh, const Expr& expr);


}


#endif
