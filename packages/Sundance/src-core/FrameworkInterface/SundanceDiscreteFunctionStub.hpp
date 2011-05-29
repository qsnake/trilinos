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

#ifndef SUNDANCE_DISCRETEFUNCTIONSTUB_H
#define SUNDANCE_DISCRETEFUNCTIONSTUB_H

#include "SundanceDefs.hpp"
#include "SundanceListExpr.hpp"
#include "SundanceDiscreteFuncDataStub.hpp"
#include "SundanceSpectralBasis.hpp"

namespace Sundance
{

/** 
 * DiscreteFunctionStub is the base class for discrete functions. 
 * Each framework will need to implement its own subclass of
 * DiscreteFunctionStub. 
 *
 * The interface is left very minimal so as to not place
 * any constraints on how a framework might specify vectors
 * and bases. When a framework needs any information about the
 * discrete function, it will have to get it by downcasting
 * to the appropriate framework-specific subclass.
 *
 * <h4> Writing a DiscreteFunctionStub subclass </h4>
 *
 * For purposes of interaction with the Sundance core, no 
 * additional methods are required.
 * However, most frameworks will require extensions to 
 * DiscreteFunctionStub that can supply the framework with information
 * on the basis and vector used by the discrete func. See the
 * demo and standard frameworks for information on how to do this.
 */
class DiscreteFunctionStub : public ListExpr
{
public:
  /** */
  DiscreteFunctionStub(const std::string& name, 
    int tensorOrder=0, int dim=1,
    const RCP<DiscreteFuncDataStub>& data=RCP<DiscreteFuncDataStub>(),
    int listIndex=0);
  /** */
  DiscreteFunctionStub(const Array<string>& name, 
    const Array<std::pair<int,int> >& tensorStructure,
    const RCP<DiscreteFuncDataStub>& data=RCP<DiscreteFuncDataStub>());
     
  /** */
  DiscreteFunctionStub(const std::string& name, 
    const SpectralBasis& sbasis, 
    int tensorOrder=0, int dim=1,
    const RCP<DiscreteFuncDataStub>& data=RCP<DiscreteFuncDataStub>(),
    int listIndex=0);
     
  /** */
  DiscreteFunctionStub(const Array<string>& name, 
    const SpectralBasis& sbasis, 
    const Array<std::pair<int,int> >& tensorStructure,
    const RCP<DiscreteFuncDataStub>& data=RCP<DiscreteFuncDataStub>());

  /** virtual destructor */
  virtual ~DiscreteFunctionStub() {;}

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}

protected:
  /** */
  void initTensor(const std::string& name, 
    int tensorOrder, int dim,
    const RCP<DiscreteFuncDataStub>& data,
    int listIndex);
  /** */
  void initTensorSpectral(const std::string& name, 
    const SpectralBasis& sbasis, 
    int tensorOrder, int dim,
    const RCP<DiscreteFuncDataStub>& data,
    int listIndex);

  /** */
  const RCP<DiscreteFuncDataStub>& dataStub() const {return data_;}

private:
  RCP<DiscreteFuncDataStub> data_;
};
}



#endif
