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
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Array.hpp"	
#include "Teuchos_Version.hpp"
#include "Optika_GUI.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerboseObject.hpp"

int main(int argc, char* argv[])
{
 /*
  * Welcome to the Optika Package!
  *
  * This package was designed to assist in the rapid development of GUIs for existing and new 
  * projects using the Trilinos Framework. Using the ParameterList class found in the Teuchos 
  * package and the new Dependency Sheet class provided by the Optika Package, this package
  * will allow you to use ParameterLists to define a set of values you wish to obtain from
  * the user. You may then pass this ParameterList to the function getInput. This function 
  * will dynamically generate a GUI based on your ParameterList, display the GUI to the user, 
  * obtain input from the user, and then store the users input back into the ParameterList. 
  * Let's take a look at an example to see how this all works.
  *
  * Before you Start:
  * We recommend you have at least a basic understanding of what a Teuchos::RCP is. While not
  * crucial to the understanding of these examples, undestanding RCPs allow you to more
  * easily understand what is going on in the examples.
  */

  /* 
   * First we create an empty parameter list. We will use this to define
   * all of the parameters we wish to obtain from the user. This type of 
   * ParameterList is commonly known as the "Valid Parameter List".
   */
  Teuchos::RCP<Teuchos::ParameterList> My_List = Teuchos::RCP<Teuchos::ParameterList>(new Teuchos::ParameterList);

  /* 
   * Creating parameters in this list can be easily done using the set function.
   * The first argument is the name of the parameter. The second is the default
   * value for the parameter. The third is a short description of what the parameter 
   * is for.
   */
  My_List->set("Max Iters", 1550, "Determines the maximum number of iterations in the solver");
  My_List->set("Tolerance", 1e-10, "The tolerance used for the convergence check");
  
  /* 
   * Validators are useful for restricting the set of values that may be used for a
   * given parameter. For the "Solver" option, we will create a validator. Here we use a 
   * StringValidator and a tuple to specify which string values
   * are valid for the "Solver" option.
   */
  Teuchos::RCP<Optika::StringValidator> solverValidator = 
     Teuchos::RCP<Optika::StringValidator>(new Optika::StringValidator(Teuchos::tuple<std::string>("GMRES", "CG", "TFQMR")));
  My_List->set( "Solver", "GMRES", "The type of solver to use.", solverValidator);

  /* 
   * The Optika Package can also handle Teuchos Arrays.
   * Here we create a Teuchos::Array object of 10 doubles
   * representing an initial guess for a linear solver.
   */
  Teuchos::Array<double> doubleArray( 10, 0.0 );
  My_List->set("Initial Guess", doubleArray, "The initial guess as a RCP to an array object.");

  /* 
   * We can also create a hieiarchy of parameters by using sublists. Here we create a sublist
   * called Prec_List. Prec_List will be contained within My_List.
   */
  Teuchos::ParameterList&
    Prec_List = My_List->sublist("Preconditioner",false,"Sublist that defines the preconditioner.");

  /*
   * Now this Prec_List can be filled with other parameters:
   */
  Prec_List.set("Type", "ILU", "The tpye of preconditioner to use");
  Prec_List.set("Drop Tolerance", 1e-3
                ,"The tolerance below which entries from the\n""factorization are left out of the factors.");

  /*
   * The getInput function starts up an Optika GUI and lets the user start to input parameter values. When the user
   * has completed their data entry, the function will finish right after all of the input values are stored in My_List.
   */
  Optika::getInput(My_List);

  /*
   * Here we can print out what the user entered in nice XML format.
   */
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::writeParameterListToXmlOStream(*My_List, *out);
  

  /*
   * A Few Final Notes
   *
   * -After calling the getInput function, any parameter in My_List has the potential to have been modified.
   *  That said, no new parameters or ParameterLists will have been added and none will have been removed.
   *
   * -The GUI can only handle certain types of parameters. They are:
   *	int
   * 	short
   * 	double
   * 	float
   * 	bool
   * 	std::string
   * 	Teuchos::Array<int>
   * 	Teuchos::Array<short>
   * 	Teuchos::Array<double>
   * 	Teuchos::Array<float>
   * 	Teuchos::Array<string>
   * If you give it a ParameterList containing a parameter that is not of one of the types specified above,
   * the parameter will still be displayed in the GUI. However, the user will not be able to modify it's value.
   *
   * That's it for now. Be sure check out the other examples to see some of the more advanced 
   * features of the Optika package. If you have any suggestions or feature requests, please send them to
   * klnusbaum@gmail.com.
   */
  return 0;
}

