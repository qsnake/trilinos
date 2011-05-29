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
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Optika_GUI.hpp"
#include "Optika_SpecificParameterEntryValidators.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Optika_StandardDependencies.hpp"
int main(){
 /*
  * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              ATTENTION              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  * !!!!   PLEASE VIEW THE BASIC EXAMPLE FIRST BEFORE READING THIS EXAMPLE. IT PROVIDES FUNDAMENTAL    !!!! 
  * !!!!   KNOWLEDGE THAT WILL BE VERY HELPFUL IN UNDERSTANDING THIS EXAMPLE.                          !!!!
  * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  */ 
  

 /*
  * Every supported data type (with the exception of bool and Teuchos::Array<bool>) has 
  * at least one validator. Lets take a look at them.
  */
 
 /*
  * The validator for number types has been templated. The validator has two constuctors.
  * The first takes three arguments:
  * 1. Minimum allowed value (inclusive).
  * 2. Maximum allowed value (inclusive).
  * 3. The step. This is the how much the value of the parameter should be changed
  * when it is told to increase or decrease in the GUI. Play around with this value a bit
  * to get a good idea of what it does. It really only has meaning when used with the GUI.
  *
  * The second constructor takes no arguments. It mearly enforces that a particular
  * parameter is of a certain number type. If you use this second constructor, no
  * minimum or maximum will be set. If you would like your validator to have a minimum 
  * and no maximum you may call the setMin function after using this constructor.
  * The same can be done with the setMax if you wish to have a maximum and no minimum.
  */

 /*
  * First we create an empty parameter list. We will use this to define
  * all of the parameters we wish to obtain from the user.
  */
  Teuchos::RCP<Teuchos::ParameterList> My_List = Teuchos::RCP<Teuchos::ParameterList>(new Teuchos::ParameterList);
  

 /*
  * Here we create a validator for an int parameter. It's minimum value is 0, and it's maximum value is 10.
  * The step is 1 (the default value for the argument).
  */
  Teuchos::RCP<Optika::EnhancedNumberValidator<int> > intVali = 
  	Teuchos::rcp(new Optika::EnhancedNumberValidator<int>(0,10));

  /*
   * We then create an int parameter and use intVali as the
   * validator argument.
   */
  My_List->set("Int", 8, "Int tester", intVali);

  /*
   * Here we create an int validator with a minimum of 0, a maximum of 100
   * and a step value of 10. Try running the example program and press the up
   * and down buttons that appear in the edit box of this parameter. That's
   * the best way to explain what the step parameter specifies.
   */
  Teuchos::RCP<Optika::EnhancedNumberValidator<int> > intStepVali =
     Teuchos::rcp(new Optika::EnhancedNumberValidator<int>(0,100,10));
  My_List->set("Step int", 10, "Step int tester", intStepVali);

  /*
   * Now suppose we wanted to make a validator for a short parameter that only
   * had a minimum, and no maximum. First we creat it.
   */
  Teuchos::RCP<Optika::EnhancedNumberValidator<short> > shortVali = 
  	Teuchos::rcp(new Optika::EnhancedNumberValidator<short>());

  /*
   * We then call the setMin function with a value of 0.
   */
  shortVali->setMin(0);

  /*
   * We then apply the validator to a short parameter.
   */
  My_List->set("Short", (short)4, "short tester", shortVali);

  /*
   * Floats and Doubles have an extra argument that can be tacked on to their constructor,
   * the precision argument. This controls how many decimals are displayed to the user in 
   * the GUI, NOT THE ACTUALL PRECISION USED TO STORE THE VALUE! Here we set the step of the
   * double validator to 1e-6 and the precision to 6 decimal places. Try running the program
   * to see it in action.
   */
  Teuchos::RCP<Optika::EnhancedNumberValidator<double> > doubleVali = 
  	Teuchos::rcp(new Optika::EnhancedNumberValidator<double>(0,20,1e-6, 6));
  My_List->set("Double", (double)4.5, "double tester", doubleVali);

  /*
   * This validator is called a StringValidator. It takes a Teuchos tuple containg strings
   * and then makes sure what every parameter it is applied to is only ever set to one of those values.
   * Note: A Teuchos StringToIntegralParameterEntryValidator would also do the trick here, but they're
   * a little more complicated to use and kind of overkill for what we're doing here.
   */
  Teuchos::RCP<Optika::StringValidator> solverValidator = Teuchos::RCP<Optika::StringValidator>(
  	new Optika::StringValidator(Teuchos::tuple<std::string>( "GMRES", "CG", "TFQMR" )));
  My_List->set("Solver", "GMRES", "The type of solver to use.", solverValidator);
  
  /*
   * The other validator that you may use is the FileNameValidator. This makes sure that the user supplies
   * a valid filename for this parameter.
   */
  Teuchos::RCP<Optika::FileNameValidator> filnameVali = 
  	Teuchos::rcp(new Optika::FileNameValidator);
  My_List->set("filename", "", "filename tester", filnameVali);

  /*
   * Array validators may also be used as well. For arrays containing numbers, simply use the ArrayNumberValidator wrapper class.
   * The ArrayNumberValidator takes an ordinary EnhancedNumberValidator as an argument for its constructor, and then uses that
   * validator to validate each entry in the array.
   *
   * If you would like to emulate the functionality of the StringValidator for an array, use the
   * ArrayStringValidator wrapper class.
   *
   * If you would like to emulate the functionality of the FileNameValidator for an array, use the ArrayFileNameValidator
   * wrapper class.
   *
   * Examples of all these are shown below.
   */
  Teuchos::Array<int> intArray(10,0);
  Teuchos::Array<std::string> stringArray(10,"Option1");
  Teuchos::Array<std::string> filenameArray(3,"~/");

  My_List->set("IntArray", intArray, "intarray tester", 
	Teuchos::RCP<Optika::ArrayNumberValidator<int> >(new Optika::ArrayNumberValidator<int>(
	  Teuchos::RCP<Optika::EnhancedNumberValidator<int> >(
	  	new Optika::EnhancedNumberValidator<int>(0,20,5)	
	))));


  Teuchos::RCP<Optika::StringValidator> optionsValidator = Teuchos::RCP<Optika::StringValidator>(
  	new Optika::StringValidator(Teuchos::tuple<std::string>("Option1", "Option2", "Option3", "Option4" )));

  My_List->set("StringArray", stringArray, "string tester", 
  	Teuchos::RCP<Optika::ArrayStringValidator>(new Optika::ArrayStringValidator(optionsValidator))); 

  Teuchos::RCP<Optika::FileNameValidator> arrayFilnameVali = Teuchos::rcp(new Optika::FileNameValidator);
  
  My_List->set("Filename Array", filenameArray, "filename array tester",
  	Teuchos::RCP<Optika::ArrayFileNameValidator>(new Optika::ArrayFileNameValidator(arrayFilnameVali)));

  /*
   * Here we print ouf the user entered values in XML format.
   */
  Optika::getInput(My_List);
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::writeParameterListToXmlOStream(*My_List, *out);

  /*
   * Final Notes:
   *
   * Supported Validators:
   * 	-EnhancedNumberValidator<S>
   * 	-StringToIntegralValidator<S>
   * 	-ArrayNumberValidator<S>
   * 	-ArrayStringValidator
   * 	-ArrayFileNameValidator
   *
   * If a validator is used that the GUI doesn't recognized, the GUI simply won't
   * use it. So if you're using validators not supported by the Optika package
   * make sure to check that your data is valid after running the getInput
   * function.
   *
   * Note: There's a validator factory class so that you can make validators
   * all quick and fast like. You should check it out if you use a lot of validators.
   * Chances are it could make your life rediculously easier.
   *
   * Remember: You default value for a parameter should be within the valid range
   * for the validator you are using on it. If this is not the case, an error will
   * be thrown before the GUI can even start up.
   *
   * Remember: Make sure you're using the right type of validator with the right
   * type of parameter, especially when it comes to arrays. If you don't use the
   * correct validator on the parameter, an error will be thrown before the GUI
   * can even startup.
   *
   * Careful! Reusing validators on multiple parameters is perfectly legal, but
   * be careful. Things like the NumberValidatorDependency can cause validators
   * min's and max's to change during the running of the GUI. Because we're
   * using pointers this means any parameter using the same validator will have
   * it's min's and max's changed too!
   */
  return 0;
}

