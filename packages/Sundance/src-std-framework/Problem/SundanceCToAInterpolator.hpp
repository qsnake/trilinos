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

#ifndef SUNDANCE_CTOAINTERPOLATOR_H
#define SUNDANCE_CTOAINTERPOLATOR_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceExpr.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "SundanceAToCPointLocator.hpp"

namespace Sundance
{
  using namespace Sundance;
  using namespace Sundance;
  using namespace Sundance;
  using namespace Sundance;
  using namespace Sundance;
  using namespace Teuchos;
  using namespace Thyra;

  /**
   * CToAInterpolator interpolates a discrete function at particle positions
   *
   * Note: not tested in parallel.
   */
  class CToAInterpolator
  {
  public:
    /** */
    CToAInterpolator(const AToCPointLocator& locator,
                     const Expr& field);

    /** */
    void interpolate(const Teuchos::Array<double>& positions,
                     Teuchos::Array<double>& results) const ;

   

    /** */
    void interpolate(const std::vector<double>& positions,
      std::vector<double>& results) const 
      {
        Teuchos::Array<double> in(positions.size());
        for (int i=0; i<in.size(); i++) in[i] = positions[i];
        Teuchos::Array<double> out;
        interpolate(in, out);
        for (int i=0; i<out.size(); i++) results[i] = out[i];
      }

    /** */
    void updateField(const Expr& field) ;


  private:

    int dim_;
    int nFacets_;
    int rangeDim_;
    RCP<Array<double> > elemToVecValuesMap_;
    AToCPointLocator locator_;
  };
}


#endif
