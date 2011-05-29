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

#ifndef SUNDANCE_SPECTRALBASISBASE_H
#define SUNDANCE_SPECTRALBASISBASE_H

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceHandleable.hpp"
#include "SundanceMap.hpp"




namespace Sundance
{
  /** Base class for spectral bases. */
  class SpectralBasisBase : public Sundance::Handleable<SpectralBasisBase>
  {
  public:
    /** Construct a basis */
    SpectralBasisBase() {;}
    
    /** virtual dtor */
    virtual ~SpectralBasisBase() {;}

    /** Return the dimension of the Spectral Basis */
    virtual int getDim() const = 0 ;

    /** Return the order of the Spectral Basis */
    virtual int getOrder() const = 0 ;

    /** Return the maximum number of terms */
    virtual int nterms() const = 0 ;
    
    /** Return the basis element stored in the basis array index */
    virtual int getElement(int i) const = 0 ;

    /** Write a description to a std::string */
    virtual std::string toString() const = 0 ;

    /** expectation operator */
    virtual double expectation(int i, int j, int k) = 0 ; 

    /** Ordering operator */
    virtual bool lessThan(const SpectralBasisBase* other) const = 0;
  };
}

#endif
