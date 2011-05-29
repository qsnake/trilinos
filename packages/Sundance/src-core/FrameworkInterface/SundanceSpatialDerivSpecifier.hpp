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

#ifndef SUNDANCE_SPATIALDERIVSPECIFIER_HPP
#define SUNDANCE_SPATIALDERIVSPECIFIER_HPP

#include "SundanceEnumTypeField.hpp"
#include "SundanceMultiIndex.hpp"

namespace Sundance
{
using namespace Sundance;

/** */
enum SpatialDerivType {IdentitySDT, PartialSDT, NormalSDT, DivSDT};



/** 
 * This class is a compact description of type of spatial derivative
 * acting on an operative function: partial derivative, divergence, 
 * or normal derivative. 
 */
class SpatialDerivSpecifier : public EnumTypeField<SpatialDerivType>
{
public:
  /** Empty ctor creates an identity operator 
   * (zeroth order partial derivative) */
  SpatialDerivSpecifier();

  /** Create a spatial derivative */
  SpatialDerivSpecifier(const MultiIndex& mi);

  /** Create a derivative of a specified type and order. */
  SpatialDerivSpecifier(const SpatialDerivType& type, int order=0);

  /** Return the multiindex of a spatial partial derivative */
  const MultiIndex& mi() const ;

  /** Return true if I am a divergence */
  bool isDivergence() const ;

  /** Return true if I am a partial derivative in a coordinate direction */
  bool isPartial() const ;

  /** Return true if I am a normal derivative */
  bool isNormal() const ;

  /** Return true if I am an identity operator */
  bool isIdentity() const ;

  /** Return the order of differentiation in the normal direction */
  int normalDerivOrder() const ;

  /** Return the order of differentiation */
  int derivOrder() const ;

  /** Write me to a std::string */
  std::string toString() const ;

  /** Comparison operator for use in sorted containers */
  bool operator<(const SpatialDerivSpecifier& other) const ;

  /** Create a new derivative that increments my multiindex by the input
   * multiindex */
  SpatialDerivSpecifier derivWrtMultiIndex(const MultiIndex& mi) const ;

private:
  MultiIndex mi_;

  int normalDerivOrder_;
};


}


namespace std
{
/** \relates SpatialDerivSpecifier */
ostream& operator<<(std::ostream& os, const Sundance::SpatialDerivSpecifier& sds);
/** \relates SpatialDerivSpecifier */
ostream& operator<<(std::ostream& os, const Sundance::SpatialDerivType& sdt);
}


#endif
