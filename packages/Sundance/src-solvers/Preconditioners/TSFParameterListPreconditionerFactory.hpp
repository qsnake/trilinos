/*@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/


#ifndef TSF_PARAMETERLIST_PRECONDITIONERFACTORY_HPP
#define TSF_PARAMETERLIST_PRECONDITIONERFACTORY_HPP

#include "SundanceDefs.hpp"
#include "TSFPreconditionerFactoryBase.hpp"
#include "TSFILUKPreconditionerFactory.hpp"
#include "TSFLinearOperatorDecl.hpp"
#include "Teuchos_ParameterList.hpp"
#include "TSFILUFactorizableOp.hpp"
#include "TSFLinearSolverBaseDecl.hpp"
#include "TSFMLOperator.hpp"

namespace TSFExtended
{
using namespace Teuchos;

/**
 * 
 */
class ParameterListPreconditionerFactory
  : public PreconditionerFactoryBase<double>
{
public:

  /** Construct with a parameter list */
  ParameterListPreconditionerFactory(const ParameterList& params)
    : PreconditionerFactoryBase<double>(), params_(params)
    {
      const std::string& pName = params_.name();
      TEST_FOR_EXCEPTION(pName != "Preconditioner", std::runtime_error,
        "expected tag=Preconditioner in parameter list " << std::endl 
        << params_);
    }

  /** virtual dtor */
  virtual ~ParameterListPreconditionerFactory(){;}
    
  /** */
  Preconditioner<double>
  createPreconditioner(const LinearOperator<double>& A) const ;

  /* Handleable boilerplate */
  GET_RCP(PreconditionerFactoryBase<double>);
private:

  ParameterList params_;
};


}

#endif




