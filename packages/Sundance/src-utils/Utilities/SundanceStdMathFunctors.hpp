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

#ifndef SUNDANCE_STDMATHFUNCTORS_H
#define SUNDANCE_STDMATHFUNCTORS_H

#include "SundanceDefs.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceUnaryFunctor.hpp"
#ifdef _MSC_VER
# include "winmath.h"
#endif

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace Sundance
{
  using namespace Teuchos;

  /** */
  class PowerFunctor : public UnaryFunctor
  {
  public:
    /** */
    PowerFunctor(const double& p);
    
    /** Evaluate power function and deriv at an array of values */ 
    virtual void eval1(const double* const x, 
              int nx, 
              double* f, 
              double* df) const ;
    /** Evaluate power function at an array of values */ 
    virtual void eval0(const double* const x, int nx, double* f) const ;

    /** Evaluate power function and first two derivs at an array of values */
    virtual void eval2(const double* const x, 
                      int nx, 
                      double* f, 
                      double* df_dx,
                      double* d2f_dxx) const ;

    /** Evaluate power function and first three derivs at an array of values */
    virtual void eval3(const double* const x, 
                      int nx, 
                      double* f, 
                      double* df_dx,
                      double* d2f_dxx,
                      double* d3f_dxxx) const ;

  private:
    double p_;
  };


  
  

  SUNDANCE_UNARY_FUNCTOR(reciprocal, StdReciprocal, "reciprocal function", 
                         NonzeroDomain(), 1.0/x[i], -f[i]*f[i], -2.0*df[i]/x[i])

    SUNDANCE_UNARY_FUNCTOR(fabs, StdFabs, "absolute value", UnboundedDomain(), ::fabs(x[i]), ((x[i]>=0.0) ? x[i] : -x[i]), 0.0)

  SUNDANCE_UNARY_FUNCTOR(sign, StdSign, "sign function", UnboundedDomain(), 
                         ((x[i]>0.0) ? 1.0 : ( (x[i]<0.0) ? -1.0 : 0.0)), 
                         0.0, 0.0)

    SUNDANCE_UNARY_FUNCTOR3(exp, StdExp, "exponential function", UnboundedDomain(), ::exp(x[i]), f[i], f[i], f[i])

    SUNDANCE_UNARY_FUNCTOR3(log, StdLog, "logarithm", PositiveDomain(), ::log(x[i]), 1.0/x[i], -df[i]*df[i], -2.0*d2f[i]*df[i])

  SUNDANCE_UNARY_FUNCTOR(sqrt, StdSqrt, "square root", PositiveDomain(), ::sqrt(x[i]), 0.5/f[i], -0.5*df[i]/x[i])

    SUNDANCE_UNARY_FUNCTOR3(sin, StdSin, "sine function", UnboundedDomain(), ::sin(x[i]), ::cos(x[i]), -f[i], -df[i])

    SUNDANCE_UNARY_FUNCTOR3(cos, StdCos, "cosine function", UnboundedDomain(), ::cos(x[i]), -::sin(x[i]), -f[i], -df[i])

  SUNDANCE_UNARY_FUNCTOR(tan, StdTan, "tangent function", UnboundedDomain(),
                         ::tan(x[i]), 1.0 + f[i]*f[i], 2.0*f[i]*df[i])

  SUNDANCE_UNARY_FUNCTOR(asin, StdASin, "inverse sine", 
                         BoundedDomain(-1.0, 1.0),
                         ::asin(x[i]), 1.0/::sqrt(1.0-x[i]*x[i]),
                         x[i]*df[i]*df[i]*df[i])

  SUNDANCE_UNARY_FUNCTOR(acos, StdACos, "inverse cosine",
                         BoundedDomain(-1.0, 1.0), 
                         ::acos(x[i]), -1.0/::sqrt(1.0-x[i]*x[i]),
                         x[i]*df[i]*df[i]*df[i])

  SUNDANCE_UNARY_FUNCTOR(atan, StdATan, "inverse tangent", 
                         UnboundedDomain(),
                         ::atan(x[i]), 1.0/(1.0 + x[i]*x[i]),
                         -2.0*x[i]*df[i]*df[i])

  SUNDANCE_UNARY_FUNCTOR(sinh, StdSinh, "hyperbolic sine",
                         UnboundedDomain(),
                         ::sinh(x[i]), ::cosh(x[i]), f[i])

  SUNDANCE_UNARY_FUNCTOR(cosh, StdCosh, "hyperbolic cosine",
                         UnboundedDomain(),
                         ::cosh(x[i]), ::sinh(x[i]), f[i])

  SUNDANCE_UNARY_FUNCTOR(tanh, StdTanh, "hyperbolic tangent",
                         UnboundedDomain(),
                         ::tanh(x[i]), 1.0 - f[i]*f[i], -2.0*f[i]*df[i])

  SUNDANCE_UNARY_FUNCTOR(asinh, StdASinh, "inverse hyperbolic sine",
                         UnboundedDomain(),
                         ::asinh(x[i]), 1.0/::sqrt(1.0 + x[i]*x[i]),
                         -x[i]*df[i]*df[i]*df[i])

  SUNDANCE_UNARY_FUNCTOR(acosh, StdACosh, "inverse hyperbolic cosine",
                         LowerBoundedDomain(1.0),
                         ::acosh(x[i]), 1.0/::sqrt(x[i]*x[i]-1.0),
                         -x[i]*df[i]*df[i]*df[i])

  SUNDANCE_UNARY_FUNCTOR(atanh, StdATanh, "inverse hyperbolic tangent",
                         BoundedDomain(-1.0, 1.0), 
                         ::atanh(x[i]), 1.0/(1.0 - x[i]*x[i]),
                         2.0*x[i]*df[i]*df[i])


}

                  
#endif  /* DOXYGEN_DEVELOPER_ONLY */  

#endif
