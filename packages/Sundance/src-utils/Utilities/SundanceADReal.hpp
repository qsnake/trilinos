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

#ifndef SUNDANCE_ADREAL_H
#define SUNDANCE_ADREAL_H


#include "SundanceDefs.hpp"
#include "SundancePoint.hpp"
#include "SundanceStdMathFunctors.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace Sundance
{
/**
 * First-order automatic differentiation of a multivariable function.
 * @author Kevin Long
 */
class ADReal
{
public:
  /** Create an ADReal with value and gradient = 0 */
  ADReal() : value_(0), gradient_(0.0, 0.0, 0.0){;}
  /** Create an ADReal object that varies linearly with
   * one coordinate direction  */
  ADReal(double value, int direction, int spatialDimension)
    : value_(value), gradient_()
    {
      if (spatialDimension==1) gradient_ = Point(0.0);
      if (spatialDimension==2) gradient_ = Point(0.0, 0.0);
      if (spatialDimension==3) gradient_ = Point(0.0, 0.0, 0.0);
      gradient_[direction] = 1.0;
    }
  /** Create a constant-valued ADReal object in a multidimensional space */
  ADReal(double value, int spatialDimension)
    : value_(value), gradient_()
    {
      if (spatialDimension==1) gradient_ = Point(0.0);
      if (spatialDimension==2) gradient_ = Point(0.0, 0.0);
      if (spatialDimension==3) gradient_ = Point(0.0, 0.0, 0.0);
    }

  /** Create an ADReal with the specified value and gradient */
  ADReal(double value, const Point& gradient)
    : value_(value), gradient_(gradient) {;}

  /** unary minus */
  ADReal operator-() const ;
  /** reflexive addition */
  ADReal& operator+=(const ADReal& other) ;
  /** reflexive subtraction */
  ADReal& operator-=(const ADReal& other) ;
  /** reflexive multiplication */
  ADReal& operator*=(const ADReal& other) ;
  /** reflexive division */
  ADReal& operator/=(const ADReal& other) ;

  /** reflexive scalar addition */
  ADReal& operator+=(const double& scalar) ;
  /** reflexive scalar subtraction */
  ADReal& operator-=(const double& scalar) ;
  /** reflexive scalar multiplication */
  ADReal& operator*=(const double& scalar) ;
  /** reflexive scalar division */
  ADReal& operator/=(const double& scalar) ;

  /** addition */
  ADReal operator+(const ADReal& other) const ;
  /** subtraction */
  ADReal operator-(const ADReal& other) const ;
  /** multiplication */
  ADReal operator*(const ADReal& other) const ;
  /** division */
  ADReal operator/(const ADReal& other) const ;

  /** scalar addition */
  ADReal operator+(const double& scalar) const ;
  /** scalar subtraction */
  ADReal operator-(const double& scalar) const ;
  /** scalar multiplication */
  ADReal operator*(const double& scalar) const ;
  /** scalar division */
  ADReal operator/(const double& scalar) const ;

  /** get the value */
  const double& value() const {return value_;}
  /** get the gradient */
  const Point& gradient() const {return gradient_;}

  void reciprocate() ;

  static double& totalFlops() {static double rtn = 0; return rtn;}


  static void addFlops(const double& flops) {totalFlops() += flops;}
private:
  double value_;
  Point gradient_;
};



/** \relates ADReal */
ADReal operator+(const double& scalar,
  const ADReal& a);
/** \relates ADReal */
ADReal operator-(const double& scalar,
  const ADReal& a);
/** \relates ADReal */
ADReal operator*(const double& scalar,
  const ADReal& a);
/** \relates ADReal */
ADReal operator/(const double& scalar,
  const ADReal& a);


inline ADReal pow(const ADReal& x, const double& y)
{
  Teuchos::RCP<UnaryFunctor> func = Teuchos::rcp(new PowerFunctor(y));
  double f;
  double df;
  double xVal = x.value();
  func->eval1(&xVal, 1, &f, &df);
  return ADReal(f, df*x.gradient());
}

#define SUNDANCE_AD_FUNCTOR(opName, functorName) \
  inline ADReal opName(const ADReal& x)\
  {\
    Teuchos::RCP<UnaryFunctor> func = Teuchos::rcp(new functorName());\
    double f;\
    double df;\
    double xVal = x.value();\
    func->eval1(&xVal, 1, &f, &df);\
    return ADReal(f, df*x.gradient());}


SUNDANCE_AD_FUNCTOR(exp, StdExp)

  SUNDANCE_AD_FUNCTOR(log, StdLog)

  SUNDANCE_AD_FUNCTOR(sqrt, StdSqrt)

  SUNDANCE_AD_FUNCTOR(sin, StdSin)

  SUNDANCE_AD_FUNCTOR(cos, StdCos)


    
  }


namespace std
{
inline ostream& operator<<(std::ostream& os, const Sundance::ADReal& x)
{
  std::cerr << "AD[" << x.value() << ", grad="<< x.gradient() << "]";
  return os;
}
}



#endif



