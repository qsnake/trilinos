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

#ifndef SUNDANCE_DISCRETEFUNCTIONDATA_H
#define SUNDANCE_DISCRETEFUNCTIONDATA_H

#include "SundanceDefs.hpp"
#include "SundanceDiscreteFuncDataStub.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "TSFVectorDecl.hpp"

namespace Sundance
{
using namespace Teuchos;

/** 
 * DiscreteFunctionData 
 */
class DiscreteFunctionData : public DiscreteFuncDataStub
{
public:
  /** */
  DiscreteFunctionData(const DiscreteSpace& space);

  /** */
  DiscreteFunctionData(const DiscreteSpace& space, 
    const TSFExtended::Vector<double>& vec);

  /** */
  DiscreteFunctionData(const DiscreteSpace& space, const double& constantValue);

  /** virtual destructor */
  virtual ~DiscreteFunctionData() {;}

  /** */
  void updateGhosts() const ;

  /** */
  void setVector(const Vector<double>& vec);

  /** */
  const Vector<double>& getVector() const {return vector_;}

  /** */
  const DiscreteSpace& discreteSpace() const {return space_;}

  /** */
  const Mesh& mesh() const {return space_.mesh();}

  /** */
  const RCP<DOFMapBase>& map() const {return space_.map();}

  /** */
  RCP<const MapStructure> getLocalValues(int cellDim, 
    const Array<int>& cellLID,
    Array<Array<double> >& localValues) const ;


  /** */
  RCP<GhostView<double> > ghostView() const 
    {updateGhosts(); return ghostView_;}

  /** */
  const BasisArray& basis() const {return space_.basis();}

  /** */
  static const DiscreteFunctionData* getData(const DiscreteFuncElement* ufe);


private:

  DiscreteSpace space_;

  Vector<double> vector_;

  mutable RCP<GhostView<double> > ghostView_;

  mutable bool ghostsAreValid_;

};
}



#endif
