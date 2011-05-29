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

#ifndef SUNDANCE_FUNCTIONIDENTIFIER_H
#define SUNDANCE_FUNCTIONIDENTIFIER_H

#include "SundanceAlgebraSpecifier.hpp"

namespace Sundance
{

/** 
 * FunctionIdentifier provides a means for distinguishing between different
 * functions and different vector components of the same function. 
 * Functions discretized with vector bases will shared a common dofID,
 * because their vector components are not independent. Functions 
 * discretized componentwise will have different IDs for each component.
 */
class FunctionIdentifier
{
public:
  /** */
  FunctionIdentifier();
  /** ctor */
  FunctionIdentifier(const AlgebraSpecifier& algSpec);
  /** ctor */
  FunctionIdentifier(const FunctionIdentifier* parent,
    const AlgebraSpecifier& componentAlgSpec);

  /** */
  std::string toString() const ;

  /** Return the ID number to be used when assigning DOFs 
      for this function */
  int dofID() const {return dofID_;}

  /** If this FID corresponds to a vector component, return the 
   *  index of the coordinate direction */
  int componentIndex() const ;

  /** Return a specification of the type of object represented, i.e.,
   * a component in a coord direction, a normal component, or a whole
   * vector. */
  const AlgebraSpecifier& algSpec() const 
    {return algSpec_;}

  /** Create a new FID representing a component of "this" vector function. */
  FunctionIdentifier createComponent(int index) const ;

  /** Create a new FID representing the normal 
   * component of "this" vector function. */
  FunctionIdentifier createNormal() const ;

  /** Comparison operator for storage in sets and maps */
  bool operator<(const FunctionIdentifier& other) const ;

  /** Equality test */
  bool operator==(const FunctionIdentifier& other) const 
    {return !(*this!=other);} 

  /** Inequality test */
  bool operator!=(const FunctionIdentifier& other) const 
    {return *this < other || other < *this;}

  /** Return true if I am a vector */
  bool isVector() const {return algSpec().isVector();}

  /** Return true if I am a coordinate component */
  bool isCoordinateComponent() const {return algSpec().isCoordinateComponent();}

  /** Return true if I am a normal component */
  bool isNormalComponent() const {return algSpec().isNormal();}

  /** Return true if I am a scalar */
  bool isScalar() const {return algSpec().isScalar();}

private:

  /** Generate a unique ID */
  static int nextID() {static int id=0; id++; return id;}

  int dofID_;

  AlgebraSpecifier algSpec_;
};

/** \relates FunctionIdentifier */
FunctionIdentifier makeFuncID(int tensorOrder);

}


namespace std
{
/** \relates FunctionIdentifier */
ostream& operator<<(std::ostream& os, const Sundance::FunctionIdentifier& fid);
}

#endif
