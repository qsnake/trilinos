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

#ifndef SUNDANCE_ATOCDENSITYSAMPLER_H
#define SUNDANCE_ATOCDENSITYSAMPLER_H

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
   * AToCDensitySampler samples a distribution of particles to compute a
   * density function on a discrete space. 
   *
   * Note: not tested in parallel.
   */
  class AToCDensitySampler
  {
  public:
    /** */
    AToCDensitySampler(const AToCPointLocator& locator,
                       const VectorType<double>& vecType);
    /** */
    AToCDensitySampler(const AToCPointLocator& locator,
                       const std::vector<double>& origin,
                       const std::vector<double>& rotationalAxis,
                       const VectorType<double>& vecType);

    /** */
    Expr sample(const std::vector<double>& positions,
                const double& particleWeight) const ;

    /** */
    Expr resetCounts() const ;

    /** */
    void addToCounts(const std::vector<double>& positions,
                     const double& particleWeight,
                     Expr density) const ;


  private:
    void init();
    Point vec2point(const std::vector<double>& x) const ;
    Point normPoint(const Point& x) const ;

    DiscreteSpace discSpace_;
    int dim_;
    Mesh mesh_;
    RCP<Array<int> > elemToVecIndexMap_;
    Expr elemWeights_;
    Vector<double> elemWeightVec_;
    AToCPointLocator locator_;
    bool isAxisymmetric_;
    Point origin_;
    Point axis_;
  };
}


#endif
