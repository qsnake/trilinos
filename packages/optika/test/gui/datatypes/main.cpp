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
#include "Teuchos_FancyOStream.hpp"
#include "Optika_GUI.hpp"

int main(){
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

  Teuchos::RCP<Teuchos::ParameterList> My_List2 = Teuchos::RCP<Teuchos::ParameterList>(new Teuchos::ParameterList);
  Teuchos::ParameterList&
    validatorList = My_List2->sublist("Validator List", false, "Validator testing\nWorking June 27th 2009");
  Teuchos::RCP<Optika::FileNameValidator> filnameVali = 
  	Teuchos::RCP<Optika::FileNameValidator>(new Optika::FileNameValidator);
  validatorList.set("filename", "", "filename tester", filnameVali);
  Teuchos::RCP<Optika::EnhancedNumberValidator<int> > intVali = 
  	Teuchos::rcp(new Optika::EnhancedNumberValidator<int>(0,10,2));
  validatorList.set("Int", 8, "Int tester", intVali);
  Teuchos::RCP<Optika::EnhancedNumberValidator<short> > shortVali = 
  	Teuchos::rcp(new Optika::EnhancedNumberValidator<short>(0,10,4));
  validatorList.set("Short", (short)4, "short tester", shortVali);
  Teuchos::RCP<Optika::EnhancedNumberValidator<float> > floatVali = 
  	Teuchos::rcp(new Optika::EnhancedNumberValidator<float>(0,20,1e-2, 6));
  validatorList.set("Float", (float)4.5, "float tester", floatVali);
  Teuchos::RCP<Optika::EnhancedNumberValidator<double> > doubleVali = 
  	Teuchos::rcp(new Optika::EnhancedNumberValidator<double>(0,20,1e-2, 6));
  validatorList.set("Double", (double)4.5, "double tester", doubleVali);
  Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
    solverValidator2 = Teuchos::rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<int>(
        Teuchos::tuple<std::string>( "GMRES", "CG", "TFQMR" )
        ,"Solver"
        )
      );
  validatorList.set(
    "Solver"
    ,"GMRES" // This will be validated by solverValidator right here!
    ,"The type of solver to use."
    ,solverValidator2
    );
  Teuchos::Array<std::string> validValues;
  validValues.append("value1");
  validValues.append("value2");
  validValues.append("value3");
  Teuchos::RCP<Optika::StringValidator> stringVali2 = Teuchos::RCP<Optika::StringValidator>(new Optika::StringValidator(validValues));
  validatorList.set("Easy String", "value1", "easy string validator tester", stringVali2);

  Teuchos::ParameterList&
    NoValiList = My_List2->sublist("No validator list",false,"sublist containing data types without validators on them for checking default behavior.");
  NoValiList.set("Int1", 8, "Int tester");
  NoValiList.set("Short1", (short)4, "short tester");
  NoValiList.set("Float1", (float)4.5, "float tester");
  NoValiList.set("Double1", (double)4.5, "double tester");
  NoValiList.set("Bool1", true);
  NoValiList.set("Bool", true);
  NoValiList.set("Free String", "fee");
  
  //Arrays
  Teuchos::RCP<Optika::StringValidator> easyStringValidator = Teuchos::RCP<Optika::StringValidator>(new Optika::StringValidator(Teuchos::tuple<std::string>("value1", "value2", "value3")));
  Teuchos::Array<int> intArray(10,0);
  Teuchos::Array<short> shortArray(10,3);
  Teuchos::Array<float> floatArray(10,4.4);
  Teuchos::Array<double> doubleArray(10, 5.5);
  Teuchos::Array<std::string> stringArray(10,"CG");
  Teuchos::Array<std::string> easyStringArray(10, "value1");
  Teuchos::Array<std::string> freestringArray(10,"Blah");
  Teuchos::Array<std::string> filenameArray(3,"/net/home/f07/klnusbau/blah.txt");
  Teuchos::ParameterList&
  	ArrayList = My_List2->sublist("Arrays", false, "sublist containing arrays.");
  ArrayList.set("IntArray", intArray, "intarray tester", Teuchos::RCP<Optika::ArrayNumberValidator<int> >(new Optika::ArrayNumberValidator<int>(intVali)));
  ArrayList.set("ShortArray", shortArray, "shortarray tester", Teuchos::RCP<Optika::ArrayNumberValidator<short> >(new Optika::ArrayNumberValidator<short>(shortVali)));
  ArrayList.set("DoubleArray", doubleArray, "doublearray tester", Teuchos::RCP<Optika::ArrayNumberValidator<double> >(new Optika::ArrayNumberValidator<double>(doubleVali)));
  ArrayList.set("FloatArray", floatArray, "floatarray tester", Teuchos::RCP<Optika::ArrayNumberValidator<float> >(new Optika::ArrayNumberValidator<float>(floatVali)));
  ArrayList.set("StringArray", stringArray, "string tester", 
  Teuchos::RCP<Optika::ArrayStringValidator>(new Optika::ArrayStringValidator(solverValidator2))); 
  ArrayList.set("EasyStringArray", easyStringArray, "testing the easy validator", Teuchos::RCP<Optika::ArrayStringValidator>(new Optika::ArrayStringValidator(easyStringValidator)));
  ArrayList.set("FreeStringArray", freestringArray, "free string array tester");
  ArrayList.set("Filename Array", filenameArray, "filename array tester",
  	Teuchos::RCP<Optika::ArrayFileNameValidator>(new Optika::ArrayFileNameValidator(filnameVali)));

  Optika::getInput(My_List2);
  Teuchos::writeParameterListToXmlOStream(*My_List2, *out);
  My_List2->print(std::cout,Teuchos::ParameterList::PrintOptions().showDoc(true).indent(2).showTypes(true));

  return 0;
}
