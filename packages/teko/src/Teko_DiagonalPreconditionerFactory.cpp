/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
//  
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//  
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//  
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//  
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//  
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission. 
//  
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
// 
// ***********************************************************************
// 
// @HEADER

*/

#include "Teko_DiagonalPreconditionerFactory.hpp"
#include "Teko_DiagonalPreconditionerOp.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MpiComm.h"
#include "EpetraExt_PointToBlockDiagPermute.h"

using Teuchos::rcp;
using Teuchos::RCP;

namespace Teko {

DiagonalPrecondState::DiagonalPrecondState(){}

/*****************************************************/


DiagonalPreconditionerFactory::DiagonalPreconditionerFactory(){
}


RCP<PreconditionerState>  DiagonalPreconditionerFactory::buildPreconditionerState() const{
   DiagonalPrecondState*  mystate = new DiagonalPrecondState(); 
   return rcp(mystate);
}


LinearOp DiagonalPreconditionerFactory::buildPreconditionerOperator(LinearOp & lo,PreconditionerState & state) const{

  // Sanity check the state
  DiagonalPrecondState *MyState = dynamic_cast<DiagonalPrecondState *> (&state);
  TEUCHOS_ASSERT(MyState != 0);

  // Get the underlying Epetra_CrsMatrix, if we have one
  RCP<const Epetra_Operator> eo=Thyra::get_Epetra_Operator(*lo);
  TEUCHOS_ASSERT(eo!=Teuchos::null);
  const Epetra_CrsMatrix *MAT=dynamic_cast<const Epetra_CrsMatrix*>(&*eo);
  TEUCHOS_ASSERT(MAT);

  // Create a new EpetraExt_PointToBlockDiagPermute for the state object, if we don't have one
  EpetraExt_PointToBlockDiagPermute *BDP;
  if(MyState->BDP_==Teuchos::null){
    BDP=new EpetraExt_PointToBlockDiagPermute(*MAT);
    BDP->SetParameters(List_);
    BDP->Compute();
    MyState->BDP_=rcp(BDP);
  }

  // Build the LinearOp object  (NTS: swapping the range and domain)
  DiagonalPreconditionerOp *MyOp=new DiagonalPreconditionerOp(MyState->BDP_,lo->domain(),lo->range());
  return rcp(MyOp);
}




void DiagonalPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList & pl){
  List_=pl;

  // Reset default to invert mode if the user hasn't specified something else
  Teuchos::ParameterList & SubList=List_.sublist("blockdiagmatrix: list");	    
  SubList.set("apply mode",SubList.get("apply mode","invert"));
}


} // end namespace Teko

