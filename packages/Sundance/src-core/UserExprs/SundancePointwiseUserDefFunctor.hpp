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

#ifndef SUNDANCE_POINTWISEUSERDEFFUNCTOR_H
#define SUNDANCE_POINTWISEUSERDEFFUNCTOR_H

#include "SundanceDefs.hpp"
#include "SundanceMap.hpp"
#include "SundanceUserDefFunctor.hpp"
#include "SundanceMultiSet.hpp"
#include "SundanceExceptions.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"



namespace Sundance
{
  using namespace Sundance;
  using namespace Teuchos;

  
  

    /**
     * PointwiseUserDefFunctor0 is an implementation of UserDefFunctor for which
     * the user writes code to evaluate the function at a single quadrature point.
     * Looping over quadrature points is done by the this class.
     */
  class PointwiseUserDefFunctor0 : public UserDefFunctor
    {
    public:
      /** ctor */
      PointwiseUserDefFunctor0(const std::string& name, int domainDim, int rangeDim) ;

      /** */
      virtual ~PointwiseUserDefFunctor0(){;}

      /** */
      void evaluationCallback(int nPoints, int maxDiffOrder,
                              const double** in,
                              double** out) const ;

      /** */
      virtual void eval0(const double* in, double* out) const = 0 ;

      /** */
      virtual int maxOrder() const {return 0;}

    private:
    };


    /**
     * PointwiseUserDefFunctor1 is an implementation of UserDefFunctor for which
     * the user writes code to evaluate the function at a single quadrature point.
     * Looping over quadrature points is done by the this class.
     */
  class PointwiseUserDefFunctor1 : public PointwiseUserDefFunctor0
    {
    public:
      /** ctor */
      PointwiseUserDefFunctor1(const std::string& name, int domainDim, int rangeDim) ;

      /** */
      virtual ~PointwiseUserDefFunctor1(){;}

      /** */
      void evaluationCallback(int nPoints, int maxDiffOrder,
                              const double** in,
                              double** out) const ;

      /** */
      virtual void eval0(const double* in, double* out) const ;

      /** */
      virtual void eval1(const double* in, double* outVals, double* outDerivs) const = 0 ;

      /** */
      virtual int maxOrder() const {return 1;}

    private:
    };


    /**
     * PointwiseUserDefFunctor2 is an implementation of UserDefFunctor for which
     * the user writes code to evaluate the function at a single quadrature point.
     * Looping over quadrature points is done by the this class.
     */
  class PointwiseUserDefFunctor2 : public PointwiseUserDefFunctor1
    {
    public:
      /** ctor */
      PointwiseUserDefFunctor2(const std::string& name, int domainDim, int rangeDim) ;

      /** */
      virtual ~PointwiseUserDefFunctor2(){;}

      /** */
      void evaluationCallback(int nPoints, int maxDiffOrder,
                              const double** in,
                              double** out) const ;

      /** */
      virtual void eval0(const double* in, double* out) const ;

      /** */
      virtual void eval1(const double* in, double* outVals, double* outDerivs) const ;

      virtual void eval2(const double* in, double* outVals, double* outDerivs,
                         double* outDerivs2) const = 0 ;

      /** */
      virtual int maxOrder() const {return 2;}

    private:
    };


}


#endif
