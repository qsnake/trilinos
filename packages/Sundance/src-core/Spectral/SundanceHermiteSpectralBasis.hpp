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

#ifndef SUNDANCE_HERMITESPECTRALBASIS_H
#define SUNDANCE_HERMITESPECTRALBASIS_H

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceMap.hpp"
#include "SundanceSpectralBasisBase.hpp"

#include "cijk.h"
#include "chaos.h"



namespace Sundance
{
using Teuchos::RCP;
using Teuchos::Array;
  /** Multidimensional Hermite spectral basis.
   * See, e.g., Ghanem and Spanos.
   *
   * \author George Saad
   */
  class HermiteSpectralBasis : public SpectralBasisBase
  {
  private:
    Array<int> basis_;
    int dim_;
    int order_;
    int maxterms_;
    RCP<cijk> cijk_;
  public:
    /** Construct a full order basis */
    HermiteSpectralBasis(int dim, int order);
    
    /** Construct a truncated basis */
    HermiteSpectralBasis(int dim, int order, int nterms); 
    
    /** Construct a basis using the specified subset of elements */
    HermiteSpectralBasis(int dim, int order, const Array<int>& basisarray); 
    

    /** Return the dim of the Spectral Basis */
    int getDim() const;

    /** Return the order of the Spectral Basis */
    int getOrder() const;

    /** Return the maximum number of terms */
    int nterms() const ;
    
    /** Return the basis element stored in the basis array index */
    int getElement(int i) const;
    
    /** expectation operator */
    double expectation(int i, int j, int k); 

    /** Write to a std::string */
    std::string toString() const ;

    /* */
    GET_RCP(SpectralBasisBase);

    /** Ordering operator */
    virtual bool lessThan(const SpectralBasisBase* other) const ;
  };
}

#endif
