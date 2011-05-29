// @HEADER
// ***********************************************************************
// 
//         Optika: A Tool For Developing Parameter Obtaining GUIs
//                Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, with Sandia Corporation, the 
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Kurtis Nusbaum (klnusbaum@gmail.com) 
// 
// ***********************************************************************
// @HEADER
#include "Optika_GUI.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_FancyOStream.hpp"
int main(){
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

  //Basic Test
  Teuchos::RCP<Teuchos::ParameterList> My_List = Teuchos::RCP<Teuchos::ParameterList>(new Teuchos::ParameterList);

  double *pointer = 0;
  My_List->set("Double pointer", pointer);
  My_List->set("Max Iters", 1550, "Determines the maximum number of iterations in the solver");
  My_List->set("Tolerance", 1e-10, "The tolerance used for the convergence check");
  
  Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
    solverValidator = Teuchos::rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<int>(
        Teuchos::tuple<std::string>( "GMRES", "CG", "TFQMR" )
        ,"Solver"
        )
      );
  My_List->set(
    "Solver"
    ,"GMRES" // This will be validated by solverValidator right here!
    ,"The type of solver to use."
    ,solverValidator
    );

   Teuchos::RCP<Optika::EnhancedNumberValidator<int> > awesomenessValidator = 
   Teuchos::RCP<Optika::EnhancedNumberValidator<int> >(new Optika::EnhancedNumberValidator<int>(0,10));
    My_List->set("Awesomeness", 5, "Rate the awesomeness!!!", awesomenessValidator);

  Teuchos::Array<double> testArray( 10, 0.0 );
  
  My_List->set("Initial Guess", testArray, "The initial guess as a RCP to an array object.");

  Teuchos::ParameterList&
    Prec_List = My_List->sublist("Preconditioner",false,"Sublist that defines the preconditioner.");

  Prec_List.set("Type", "ILU", "The tpye of preconditioner to use");
  Prec_List.set("Drop Tolerance", 1e-3
                ,"The tolerance below which entries from the\n""factorization are left out of the factors.");

  Optika::getInput(My_List);

  Teuchos::writeParameterListToXmlOStream(*My_List, *out);

  return 0;
}
