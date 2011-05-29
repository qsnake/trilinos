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

#include <math.h>
#include "SundanceStdMathFunctors.hpp"
#include "Teuchos_Utils.hpp"

using namespace Sundance;
using namespace Teuchos;


PowerFunctor::PowerFunctor(const double& p) 
  : UnaryFunctor("pow("+Teuchos::toString(p)+")"), p_(p)
{;}

void PowerFunctor::eval1(const double* const x, 
                        int nx, 
                        double* f, 
                        double* df) const
{
  if (p_==2)
    {
      for (int i=0; i<nx; i++) 
        {
          df[i] = 2.0*x[i];
          f[i] = x[i]*x[i];
        }
    }
  else if (p_==1)
    {
      for (int i=0; i<nx; i++) 
        {
          df[i] = 0.0;
          f[i] = 1.0;
        }
    }
  else if (p_==0)
    {
      for (int i=0; i<nx; i++) 
        {
          df[i] = 0.0;
          f[i] = 0.0;
        }
    }
  else
    {
      if (checkResults())
        {
          for (int i=0; i<nx; i++) 
            {
              double px = ::pow(x[i], p_-1);
              df[i] = p_*px;
              f[i] = x[i]*px;
              //bvbw tried to include math.h, without success
#ifdef REDDISH_PORT_PROBLEM
              TEST_FOR_EXCEPTION(fpclassify(f[i]) != FP_NORMAL 
                                 || fpclassify(df[i]) != FP_NORMAL,
                                 RuntimeError,
                                 "Non-normal floating point result detected in "
                                 "evaluation of unary functor " << name());
#endif
            }
        }
      else
        {
          for (int i=0; i<nx; i++) 
            {
              double px = ::pow(x[i], p_-1);
              df[i] = p_*px;
              f[i] = x[i]*px;
            }
        }
    }
}

void PowerFunctor::eval2(const double* const x, 
                        int nx, 
                        double* f, 
                        double* df,
                        double* d2f_dxx) const
{
  if (p_==2)
    {
      for (int i=0; i<nx; i++) 
        {
          d2f_dxx[i] = 2.0;
          df[i] = 2.0*x[i];
          f[i] = x[i]*x[i];
        }
    }
  else if (p_==1)
    {
       for (int i=0; i<nx; i++) 
        {
          d2f_dxx[i] = 0.0;
          df[i] = 1.0;
          f[i] = x[i];
        }
    }
  else if (p_==0)
    {
      for (int i=0; i<nx; i++) 
        {
          d2f_dxx[i] = 0.0;
          df[i] = 0.0;
          f[i] = 1.0;
        }
    }
  else
    {
      if (checkResults())
        {
          for (int i=0; i<nx; i++) 
            {
              double px = ::pow(x[i], p_-2);
              d2f_dxx[i] = p_ * (p_-1) * px;
              df[i] = p_*x[i]*px;
              f[i] = x[i]*x[i]*px;
#ifdef REDDISH_PORT_PROBLEM
              TEST_FOR_EXCEPTION(fpclassify(f[i]) != FP_NORMAL 
                                 || fpclassify(df[i]) != FP_NORMAL,
                                 RuntimeError,
                                 "Non-normal floating point result detected in "
                                 "evaluation of unary functor " << name());
#endif
            }
        }
      else
        {
          for (int i=0; i<nx; i++) 
            {
              double px = ::pow(x[i], p_-2);
              d2f_dxx[i] = p_ * (p_-1) * px;
              df[i] = p_*x[i]*px;
              f[i] = x[i]*x[i]*px;
            }
        }
    }
}


void PowerFunctor::eval3(const double* const x, 
                         int nx, 
                         double* f, 
                         double* df,
                         double* d2f_dxx,
                         double* d3f_dxxx) const
{
  if (p_==3)
    {
      for (int i=0; i<nx; i++) 
        {
          d3f_dxxx[i] = 6.0;
          d2f_dxx[i] = 6.0*x[i];
          df[i] = 3.0*x[i]*x[i];
          f[i] = x[i]*x[i]*x[i];
        }
    }
  else if (p_==2)
    {
      for (int i=0; i<nx; i++) 
        {
          d3f_dxxx[i] = 0.0;
          d2f_dxx[i] = 2.0;
          df[i] = 2.0*x[i];
          f[i] = x[i]*x[i];
        }
    }
  else if (p_==1)
    {
       for (int i=0; i<nx; i++) 
        {
          d3f_dxxx[i] = 0.0;
          d2f_dxx[i] = 0.0;
          df[i] = 1.0;
          f[i] = x[i];
        }
    }
  else if (p_==0)
    {
      for (int i=0; i<nx; i++) 
        {
          d3f_dxxx[i] = 0.0;
          d2f_dxx[i] = 0.0;
          df[i] = 0.0;
          f[i] = 1.0;
        }
    }
  else
    {
      if (checkResults())
        {
          for (int i=0; i<nx; i++) 
            {
              double px = ::pow(x[i], p_-3);
              d3f_dxxx[i] = p_ * (p_-1) * (p_-2) * px;
              d2f_dxx[i] = p_ * (p_-1) * x[i] * px;
              df[i] = p_*x[i]*x[i]*px;
              f[i] = x[i]*x[i]*x[i]*px;
#ifdef REDDISH_PORT_PROBLEM
              TEST_FOR_EXCEPTION(fpclassify(f[i]) != FP_NORMAL 
                                 || fpclassify(df[i]) != FP_NORMAL,
                                 RuntimeError,
                                 "Non-normal floating point result detected in "
                                 "evaluation of unary functor " << name());
#endif
            }
        }
      else
        {
          for (int i=0; i<nx; i++) 
            {
              double px = ::pow(x[i], p_-3);
              d3f_dxxx[i] = p_ * (p_-1) * (p_-2) * px;
              d2f_dxx[i] = p_ * (p_-1) * x[i] * px;
              df[i] = p_*x[i]*x[i]*px;
              f[i] = x[i]*x[i]*x[i]*px;
            }
        }
    }
}

void PowerFunctor::eval0(const double* const x, 
                        int nx, 
                        double* f) const
{
  if (checkResults())
    {
      for (int i=0; i<nx; i++) 
        {
          f[i] = ::pow(x[i], p_);
#ifdef REDDISH_PORT_PROBLEM
          TEST_FOR_EXCEPTION(fpclassify(f[i]) != FP_NORMAL, 
                             RuntimeError,
                             "Non-normal floating point result detected in "
                             "evaluation of unary functor " << name());
#endif
        }
    }
  else
    {
      for (int i=0; i<nx; i++) 
        {
          f[i] = ::pow(x[i], p_);
        }
    }
}
