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

#ifndef SUNDANCE_ADFIELD_H
#define SUNDANCE_ADFIELD_H

#include "SundanceDefs.hpp"
#include "SundanceADBasis.hpp"

namespace SundanceTesting
{
  using namespace Sundance;
  using namespace Teuchos;
  using namespace Sundance;
  using namespace Sundance;

  /** 
   *
   */
  class ADField
  {
  public:
    /** */
    ADField(){;}

    /** */
    ADField(const ADBasis& basis, const double& coeff);

    /** */
    static Point& evalPoint() {static Point x(0.123, 0.456, 0.789); return x;}


    /** */
    void setCoeff(const double& c);

    /** */
    const ADBasis& basis() const {return basis_;}

    /** */
    ADReal evaluate() const ;

    /** */
    double coeff() const {return *coeff_;}

    double value() const {return evaluate().value();}

    ADReal operator+(const ADReal& x) const ;

    ADReal operator+(const double& x) const ;

    ADReal operator+(const ADField& x) const ;



    ADReal operator-(const ADReal& x) const ;

    ADReal operator-(const double& x) const ;

    ADReal operator-(const ADField& x) const ;



    ADReal operator*(const ADReal& x) const ;

    ADReal operator*(const double& x) const ;

    ADReal operator*(const ADField& x) const ;



    ADReal operator/(const ADReal& x) const ;

    ADReal operator/(const double& x) const ;

    ADReal operator/(const ADField& x) const ;



    ADReal operator-() const ;

    ADReal reciprocal() const ;
    

    
  private:
    ADBasis basis_;
    RCP<double> coeff_;
  };


  inline ADReal operator+(const ADReal& x, const ADField& y)
  {
    return y + x;
  }

  inline ADReal operator+(const double& x, const ADField& y)
  {
    return y + x;
  }

  inline ADReal operator-(const ADReal& x, const ADField& y)
  {
    return -y + x;
  }

  inline ADReal operator-(const double& x, const ADField& y)
  {
    return -y + x;
  }

  inline ADReal operator*(const ADReal& x, const ADField& y)
  {
    return y * x;
  }

  inline ADReal operator*(const double& x, const ADField& y)
  {
    return y * x;
  } 

  inline ADReal operator/(const ADReal& x, const ADField& y)
  {
    return x * y.reciprocal();
  }

  inline ADReal operator/(const double& x, const ADField& y)
  {
    return x * y.reciprocal();
  } 
  
  inline ADReal sin(const ADField& x)
  {
    return sin(x.evaluate());
  }
  
  inline ADReal cos(const ADField& x)
  {
    return cos(x.evaluate());
  }
  
  inline ADReal pow(const ADField& x, const double& y)
  {
    return pow(x.evaluate(), y);
  }
  
  inline ADReal exp(const ADField& x)
  {
    return exp(x.evaluate());
  }
  
  inline ADReal log(const ADField& x)
  {
    return log(x.evaluate());
  }
}



#endif
