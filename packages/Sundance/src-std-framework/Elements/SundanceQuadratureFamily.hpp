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

#ifndef SUNDANCE_QUADRATUREFAMILY_H
#define SUNDANCE_QUADRATUREFAMILY_H

#include "SundanceDefs.hpp"
#include "SundanceQuadratureFamilyBase.hpp"
#include "SundanceHandle.hpp"

namespace Sundance
{

/** 
 * QuadratureFamily is a geometry-independent specification of
 * a method by which quadrature is to be carried out. For example,
 * a GaussianQuadrature family will generate Gaussian
 * quadrature points on any cell type.
 */
class QuadratureFamily : public Sundance::Handle<QuadratureFamilyStub>
{
public:
  /* */
  HANDLE_CTORS(QuadratureFamily, QuadratureFamilyStub);
  /** */
  XMLObject toXML() const ;

  /** */
  int order() const ;

  /** Returns the number of points in a rule of the given cell type 
      WARNING: this is slow.  Call it once and store the result. 
      TODO: make it pure virtual and override with queries in
      the derived classes, making them supply the information.  */
  int getNumPoints( const CellType& cellType ) const;


  /** Get the quadrature points and weights for the given cell type */
  void getPoints(const CellType& cellType, 
    Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;

  /** Get quadrature points and weights for integration on a facet of a cell */
  void getFacetPoints(const CellType& cellType, 
    int facetDim,
    int facetIndex,
    Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;

  /** Get the quadrature points and weights for the given cell type ,
   * which might be cut by a curve in the case of, see class for more info */
  void getAdaptedWeights(const CellType& cellType ,
  		         int cellDim,
	             int celLID,
	             int facetIndex ,
                 const Mesh& mesh ,
                 const ParametrizedCurve& globalCurve ,
                 Array<Point>& quadPoints ,
                 Array<double>& quadWeights ,
                 bool& isCut) const ;

private:
  /** Get quad points for a facet of a line */
  void getLineFacetQuad(int facetDim,
    int facetIndex,
    Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;
  /** Get quad points for a facet of a triangle */
  void getTriangleFacetQuad(int facetDim,
    int facetIndex,
    Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;

  /** Get quad points for a facet of a quadlateral */
  void getQuadFacetQuad(int facetDim,
    int facetIndex,
    Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;

  /** Get quad points for a facet of a tet */
  void getTetFacetQuad(int facetDim,
    int facetIndex,
    Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;

  /** Get quad points for a facet of a Brick cell */
  void getBrickFacetQuad(int facetDim,
    int facetIndex,
    Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;
};


/** \relates QuadratureFamily */
void printQuad(std::ostream& os, 
  const Array<Point>& pts, const Array<double>& wgts);


}

#endif
