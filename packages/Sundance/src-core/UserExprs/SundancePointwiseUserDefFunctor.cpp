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

#include "SundancePointwiseUserDefFunctor.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;


PointwiseUserDefFunctor0::PointwiseUserDefFunctor0(const std::string& name, 
                                                   int domainDim, 
                                                   int rangeDim)
  :  UserDefFunctor(name, domainDim, rangeDim)
{}

void PointwiseUserDefFunctor0::evaluationCallback(int nPoints, int maxDiffOrder,
                                    const double** in,
                                    double** out) const 
{
  TEST_FOR_EXCEPTION(maxDiffOrder > 0, RuntimeError,
                     "diff order = " << maxDiffOrder 
                     << " not supported for functor "
                     << name());

  static Array<double> x;
  x.resize(domainDim());
  
  static Array<double> f;
  f.resize(rangeDim());

  for (int i=0; i<nPoints; i++)
    {
      for (int j=0; j<domainDim(); j++) x[j] = in[j][i];
      const double* xp = &(x[0]);
      double* fp = &(f[0]);
      eval0(xp, fp);
      for (int j=0; j<rangeDim(); j++) out[j][i] = fp[j];
    }
}

PointwiseUserDefFunctor1::PointwiseUserDefFunctor1(const std::string& name, 
                                                   int domainDim, 
                                                   int rangeDim)
  : PointwiseUserDefFunctor0(name, domainDim, rangeDim)
{}

void PointwiseUserDefFunctor1::evaluationCallback(int nPoints, int maxDiffOrder,
                                    const double** in,
                                    double** out) const 
{
  TEST_FOR_EXCEPTION(maxDiffOrder > 1, RuntimeError,
                     "diff order = " << maxDiffOrder 
                     << " not supported for functor "
                     << name());

  static Array<double> x;
  x.resize(domainDim());
  
  static Array<double> f;
  if (maxDiffOrder==1) f.resize(rangeDim() * (1 + domainDim()) );
  else f.resize(rangeDim());
  double* fp = &(f[0]);

  for (int i=0; i<nPoints; i++)
    {
      for (int j=0; j<domainDim(); j++) x[j] = in[j][i];
      const double* xp = &(x[0]);

      if (maxDiffOrder==1) 
        {
          double* dfp = &(f[rangeDim()]);
          eval1(xp, fp, dfp);
        }
      else eval0(xp, fp);
      for (int j=0; j<f.size(); j++) out[j][i] = fp[j];
    }
}


void PointwiseUserDefFunctor1::eval0(const double* in, double* out) const 
{
  static Array<double> dummy;
  dummy.resize(domainDim() * rangeDim());

  eval1(in, out, &(dummy[0]));
}



PointwiseUserDefFunctor2::PointwiseUserDefFunctor2(const std::string& name, 
                                                   int domainDim, 
                                                   int rangeDim)
  : PointwiseUserDefFunctor1(name, domainDim, rangeDim)
{}

void PointwiseUserDefFunctor2::evaluationCallback(int nPoints, int maxDiffOrder,
                                                  const double** in,
                                                  double** out) const 
{
  TEST_FOR_EXCEPTION(maxDiffOrder > 2 || maxDiffOrder < 0, RuntimeError,
                     "diff order = " << maxDiffOrder 
                     << " not supported for functor "
                     << name());

  int nTotal = 1;
  int numFirst = domainDim();
  int numSecond = domainDim()*(domainDim()+1)/2;

  static Array<double> x;
  x.resize(domainDim());
  
  static Array<double> f;

  if (maxDiffOrder > 0) nTotal += numFirst;
  if (maxDiffOrder > 1) nTotal += numSecond;

  f.resize(rangeDim() * nTotal);

  double* fp = &(f[0]);

  for (int i=0; i<nPoints; i++)
    {
      for (int j=0; j<domainDim(); j++) x[j] = in[j][i];
      const double* xp = &(x[0]);

      if (maxDiffOrder==0)
        {
          eval0(xp, fp);
        }
      else if (maxDiffOrder==1)
        {
          double* dfp = &(f[rangeDim()]);
          eval1(xp, fp, dfp);
        }
      else if (maxDiffOrder==2)
        {
          double* dfp = &(f[rangeDim()]);
          double* d2fp = &(f[rangeDim()*(1 + domainDim())]);
          eval2(xp, fp, dfp, d2fp);
        }
      else
        {
          TEST_FOR_EXCEPT(true);
        }
      for (int j=0; j<f.size(); j++) out[j][i] = fp[j];
    }
}


void PointwiseUserDefFunctor2::eval0(const double* in, double* f) const 
{
  static Array<double> dummy1;
  static Array<double> dummy2;
  dummy1.resize(rangeDim() *  domainDim());
  dummy2.resize(rangeDim() *  domainDim()*(domainDim()+1)/2);

  eval2(in, f, &(dummy1[0]), &(dummy2[0]));
}


void PointwiseUserDefFunctor2::eval1(const double* in, double* f, double* df) const 
{
  static Array<double> dummy2;
  dummy2.resize(rangeDim() *  domainDim()*(domainDim()+1)/2);

  eval2(in, f, df, &(dummy2[0]));
}
