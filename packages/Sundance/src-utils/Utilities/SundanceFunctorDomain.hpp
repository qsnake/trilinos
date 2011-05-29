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

#ifndef SUNDANCE_FUNCTORDOMAIN_H
#define SUNDANCE_FUNCTORDOMAIN_H

#include "SundanceDefs.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace Sundance
{
  using namespace Teuchos;

  class FunctorDomain
  {
  public:
    FunctorDomain();

    virtual ~FunctorDomain(){;}

    virtual bool hasLowerBound() const {return false;}

    virtual double lowerBound() const ;

    virtual bool hasUpperBound() const {return false;}

    virtual double upperBound() const ;

    virtual bool hasExcludedPoint() const {return false;}

    virtual double excludedPoint() const ;

  };

  class UnboundedDomain : public FunctorDomain
  {
  public:
    UnboundedDomain();
  };


  class PositiveDomain : public FunctorDomain
  {
  public:
    PositiveDomain();

     bool hasLowerBound() const {return true;}

     double lowerBound() const {return 0.0;}
  };


  class BoundedDomain : public FunctorDomain
  {
  public:
    BoundedDomain(const double& lower, const double& upper);

     bool hasLowerBound() const {return true;}

     double lowerBound() const {return lower_;}

     bool hasUpperBound() const {return true;}

     double upperBound() const {return upper_;}

  private:
    double lower_;

    double upper_;
  };


  class LowerBoundedDomain : public FunctorDomain
  {
  public:
    LowerBoundedDomain(const double& lower);

     bool hasLowerBound() const {return true;}

     double lowerBound() const {return lower_;}

  private:
    double lower_;
  };

class NonzeroDomain : public FunctorDomain
  {
  public:
    NonzeroDomain();

     bool hasExcludedPoint() const {return true;}

     double excludedPoint() const {return 0.0;}
  };

}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
