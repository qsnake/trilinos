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

#ifndef SUNDANCE_SUMTRANSFORMATIONSEQUENCE_H
#define SUNDANCE_SUMTRANSFORMATIONSEQUENCE_H

#include "SundanceDefs.hpp"
#include "SundanceSumTransformation.hpp"
#include "Teuchos_Array.hpp"

namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;




/** 
 * SumTransformationSequence is a sequence of transformations
 * to be applied to a sum, producing a transformed expression. 
 */
class SumTransformationSequence : public SumTransformation,
                                  public Array<RCP<SumTransformation> >
{
public:
  /** */
  SumTransformationSequence();

  /** */
  virtual ~SumTransformationSequence(){;}

  /**
   * Test whether the transform is applicable in this case,
   * and if it is, apply it. The return value is true is the
   * transformation was applied, otherwise false. 
   * Returns by non-const reference
   * the transformed expression. 
   *
   * For SumTransformationSequence, this is implemented by
   * trying to apply all transformations in sequence. If one
   * succeeds, we exit immediately with true.
   */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    int sign, RCP<ScalarExpr>& rtn) const ;

        
};
}

#endif
