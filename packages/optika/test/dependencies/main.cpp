// @HEADER // ***********************************************************************
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
#include "Teuchos_LocalTestingHelpers.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Optika_StandardDependencies.hpp"
#include "Optika_DependencySheet.hpp"
#include "Optika_StandardConditions.hpp"

int intFuncTester(int argument){
	return argument+10;
}

int intVisualTester(int argument){
	if(argument <= 32){
		return 1;
	}
	else{
		return 0;
	}
}

double fondueTempTester(double argument){
	return argument-100.0;
}


/**
 * Test all the validator dependencies.
 */
int testValiDeps(Teuchos::FancyOStream &out){
	bool success = true;
	Teuchos::RCP<Teuchos::ParameterList> My_deplist = Teuchos::RCP<Teuchos::ParameterList>(new Teuchos::ParameterList);
	Teuchos::RCP<Optika::DependencySheet> depSheet1 = Teuchos::RCP<Optika::DependencySheet>(new Optika::DependencySheet(My_deplist));

	/*
	 * Testing StringValidatorDependency
	 */
 	Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
   	stringFoodTypeValidator = Teuchos::rcp(
		new Teuchos::StringToIntegralParameterEntryValidator<int>(
		Teuchos::tuple<std::string>( "Cheese", "Soda", "Chips" )
		,"Food Type"
		)
	);

	Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
    cheeseValidator = Teuchos::rcp(
		new Teuchos::StringToIntegralParameterEntryValidator<int>(
   			Teuchos::tuple<std::string>( "Swiss", "American", "Super Awesome Cheese" )
			,"Food Selector"
			)
	);

	Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
	sodaValidator = Teuchos::rcp(
		new Teuchos::StringToIntegralParameterEntryValidator<int>(
			Teuchos::tuple<std::string>( "Pepsi", "Coke", "Kurtis Cola", "Bad Cola" )
			,"Food Selector"
			)
		);

	Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
	chipsValidator = Teuchos::rcp(
		new Teuchos::StringToIntegralParameterEntryValidator<int>(
			Teuchos::tuple<std::string>( "Lays", "Doritos", "Kurtis Super Awesome Brand" )
			,"Food Selector"
		)
	);

	Optika::StringValidatorDependency::ValueToValidatorMap testValidatorMap1;
	testValidatorMap1["Cheese"] = cheeseValidator;
	testValidatorMap1["Soda"] = sodaValidator;
	testValidatorMap1["Chips"] = chipsValidator;

	Teuchos::ParameterList& stringValiDepList = My_deplist->sublist("String Validator Dependency", false, "String Validator Dependency testing list.\nWorking June 27th 2009");
	stringValiDepList.set("Food Selector", "Swiss", "select the food you want", cheeseValidator);
	stringValiDepList.set("Food Type", "Cheese", "String Validator Dependency Tester", stringFoodTypeValidator);

	Teuchos::RCP<Optika::StringValidatorDependency> 
	stringValiDep = Teuchos::RCP<Optika::StringValidatorDependency>(
		new Optika::StringValidatorDependency(
			"Food Type", 
			Teuchos::sublist(My_deplist,"String Validator Dependency"),
			"Food Selector", 
			Teuchos::sublist(My_deplist,"String Validator Dependency"),
			testValidatorMap1, 
			cheeseValidator
		)
	);

	depSheet1->addDependency(stringValiDep);
	
	TEST_NOTHROW(stringValiDepList.validateParameters(stringValiDepList));
	TEST_ASSERT(depSheet1->hasDependents(stringValiDepList.getEntryPtr("Food Type")));
	Optika::DependencySheet::DepSet stringValiDepSet = depSheet1->getDependenciesForParameter(stringValiDepList.getEntryPtr("Food Type"));
	TEST_ASSERT(stringValiDepSet.size() == 1);
	stringValiDepList.set("Food Type","Soda");
	stringValiDep->evaluate();
	TEST_ASSERT(stringValiDepList.getEntry("Food Selector").validator().get()==sodaValidator.get());
	TEST_THROW(stringValiDepList.validateParameters(stringValiDepList), Teuchos::Exceptions::InvalidParameterValue);
	stringValiDepList.set("Food Selector", "Pepsi");
	TEST_NOTHROW(stringValiDepList.validateParameters(stringValiDepList));
 

	/*
	 * Tesing some different aspects of the StringValidatorDependency
	 */
	Teuchos::ParameterList& 
	stringValiDepList2 = My_deplist->sublist(
		"String Validator Dependency (other validators)",
		false,
		"String Validator Dependency testing list for EnhancedNumber Validators.\nWorking June 27th 2009"
	);

	Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
	stringRangeValidator = Teuchos::rcp(
		new Teuchos::StringToIntegralParameterEntryValidator<int>(
		Teuchos::tuple<std::string>( "1-10", "10-33", "50-60" ),
		"Range selector"
		)
	);

	Teuchos::RCP<Optika::EnhancedNumberValidator<int> > range110Vali = 
	Teuchos::rcp(new Optika::EnhancedNumberValidator<int>(1,10));
	Teuchos::RCP<Optika::EnhancedNumberValidator<int> > range1033Vali = 
	Teuchos::rcp(new Optika::EnhancedNumberValidator<int>(10,33));
	Teuchos::RCP<Optika::EnhancedNumberValidator<int> > range5060Vali = 
	Teuchos::rcp(new Optika::EnhancedNumberValidator<int>(50,60));

	stringValiDepList2.set("Range selector", "1-10", "selects the range to validate", stringRangeValidator);

	Optika::StringValidatorDependency::ValueToValidatorMap rangeValidatorMap1;
	rangeValidatorMap1["1-10"] = range110Vali;
	rangeValidatorMap1["10-33"] = range1033Vali;
	rangeValidatorMap1["50-60"] = range5060Vali;
	stringValiDepList2.set("RangeValue", 3, "the value of the range", range110Vali);

	Teuchos::RCP<Optika::StringValidatorDependency> 
	stringValiDep2 = Teuchos::RCP<Optika::StringValidatorDependency>(
		new Optika::StringValidatorDependency(
			"Range selector", 
			Teuchos::sublist(My_deplist,"String Validator Dependency (other validators)"),
			"RangeValue", 
			Teuchos::sublist(My_deplist,"String Validator Dependency (other validators)"),
			rangeValidatorMap1, 
			range110Vali
		)
	);

	depSheet1->addDependency(stringValiDep2);

	TEST_NOTHROW(stringValiDepList2.validateParameters(stringValiDepList2));
	TEST_ASSERT(depSheet1->hasDependents(stringValiDepList2.getEntryPtr("Range selector")));
	Optika::DependencySheet::DepSet stringValiDepSet2 = depSheet1->getDependenciesForParameter(stringValiDepList2.getEntryPtr("Range selector"));
	TEST_ASSERT(stringValiDepSet2.size() == 1);
	stringValiDepList2.set("Range selector","50-60");
	stringValiDep2->evaluate();
	TEST_ASSERT(stringValiDepList2.getEntry("RangeValue").validator().get() == range5060Vali.get());
	TEST_THROW(stringValiDepList2.validateParameters(stringValiDepList2), Teuchos::Exceptions::InvalidParameterValue);
	stringValiDepList2.set("RangeValue", 55);
	TEST_NOTHROW(stringValiDepList2.validateParameters(stringValiDepList2));

	/*
	 * Testing the BoolValidatorDependency.
	 */
	Teuchos::ParameterList&
	boolValidatorDepList = My_deplist->sublist(
		"Bool Validator Dependency List", 
		false,
		"Bool Validator Dependency testing list.\nWorking June 27th 2009"
	);

	boolValidatorDepList.set("Use Validator?", true, "truns the validator on and off");
	Teuchos::RCP<Optika::EnhancedNumberValidator<int> > basicVali = Teuchos::rcp(new Optika::EnhancedNumberValidator<int>(1,10));
	Teuchos::RCP<Optika::EnhancedNumberValidator<int> > basicVali2 = Teuchos::rcp(new Optika::EnhancedNumberValidator<int>());
	boolValidatorDepList.set("do I have a validator?", 4, "does it have a validator?", basicVali);

	Teuchos::RCP<Optika::BoolValidatorDependency> 
	boolValiDep = Teuchos::RCP<Optika::BoolValidatorDependency>(
		new Optika::BoolValidatorDependency(
			"Use Validator?", 
			Teuchos::sublist(My_deplist,"Bool Validator Dependency List"),
			"do I have a validator?", 
			Teuchos::sublist(My_deplist,"Bool Validator Dependency List"),
			basicVali, 
			basicVali2
		)
	);

	depSheet1->addDependency(boolValiDep);

	TEST_ASSERT(depSheet1->hasDependents(boolValidatorDepList.getEntryPtr("Use Validator?")));
	TEST_ASSERT(boolValidatorDepList.getEntry("do I have a validator?").validator().get() == basicVali.get());
	TEST_NOTHROW(boolValidatorDepList.validateParameters(boolValidatorDepList));
	Optika::DependencySheet::DepSet boolValiDepSet = depSheet1->getDependenciesForParameter(boolValidatorDepList.getEntryPtr("Use Validator?"));
	TEST_ASSERT(boolValiDepSet.size() == 1);
	boolValidatorDepList.set("Use Validator?",false);
	boolValiDep->evaluate();
	TEST_ASSERT(boolValidatorDepList.getEntry("do I have a validator?").validator().get() == basicVali2.get());


	/*
	 * Testing the RangeValidatorDependency
	 */
	Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
	lowTempCheeseValidator = Teuchos::rcp(
		new Teuchos::StringToIntegralParameterEntryValidator<int>(
			Teuchos::tuple<std::string>( "PepperJack", "Swiss", "American" ),
			"Cheese to Fondue"
		)
	);

	Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
	highTempCheeseValidator = Teuchos::rcp(
		new Teuchos::StringToIntegralParameterEntryValidator<int>(
			Teuchos::tuple<std::string>( "Munster", "Provalone", "Kurtis Super Awesome Cheese"),
			"Cheese to Fondue"
		)
	);

	Teuchos::ParameterList& 
	rangeValidatorDepList = My_deplist->sublist(
		"Range Validator Dependency List",
		false,
		"Range Validator Dependency testing list.\nWorking June 27th 2009"
	);
	rangeValidatorDepList.set("Temperature",101.0, "The temperature of the fondue");
	rangeValidatorDepList.set("Cheese to Fondue", "Swiss", "The cheese we'll be using in our fondue pot.", lowTempCheeseValidator);
	Optika::RangeValidatorDependency<double>::RangeToValidatorMap tempranges;
	tempranges[std::pair<double,double>(100,200)] = lowTempCheeseValidator;
	tempranges[std::pair<double,double>(200,300)] = highTempCheeseValidator;
	Teuchos::RCP<Optika::RangeValidatorDependency<double> > 
	cheeseTempDep = Teuchos::RCP<Optika::RangeValidatorDependency<double> >(
		new Optika::RangeValidatorDependency<double>(
			"Temperature", 
			Teuchos::sublist(My_deplist,"Range Validator Dependency List"),
			"Cheese to Fondue", 
			Teuchos::sublist(My_deplist,"Range Validator Dependency List"),
			tempranges, 
			lowTempCheeseValidator
		)
	);
	depSheet1->addDependency(cheeseTempDep);

	TEST_ASSERT(depSheet1->hasDependents(rangeValidatorDepList.getEntryPtr("Temperature")));
	Optika::DependencySheet::DepSet rangeValiDepSet = depSheet1->getDependenciesForParameter(rangeValidatorDepList.getEntryPtr("Temperature"));
	TEST_ASSERT(rangeValiDepSet.size() == 1);
	rangeValidatorDepList.set("Temperature",250.0);
	cheeseTempDep->evaluate();
	TEST_ASSERT(rangeValidatorDepList.getEntry("Cheese to Fondue").validator().get() == highTempCheeseValidator.get());
	TEST_THROW(rangeValidatorDepList.validateParameters(rangeValidatorDepList), Teuchos::Exceptions::InvalidParameterValue);
	rangeValidatorDepList.set("Cheese to Fondue", "Provalone");
	TEST_NOTHROW(rangeValidatorDepList.validateParameters(rangeValidatorDepList));

	return (success ? 0:1);
}
  
/**
 * Testing all the visual dependencies
 */
int testVisualDeps(Teuchos::FancyOStream &out){
	bool success = true;
	Teuchos::RCP<Teuchos::ParameterList> My_deplist = Teuchos::RCP<Teuchos::ParameterList>(new Teuchos::ParameterList);
	Teuchos::RCP<Optika::DependencySheet> depSheet1 = Teuchos::RCP<Optika::DependencySheet>(new Optika::DependencySheet(My_deplist));

	/*
	 * Testing the NumberVisualDependency
	 */
	Teuchos::ParameterList&
	doubleVisualDepList = My_deplist->sublist(
		"NumberVisual Dependency List (double)", 
		false, 
		"Number visual Dependency testing list.\nWorking June 27th 2009"
	);
		
	doubleVisualDepList.set("Temperature",101.0, "The temperature of the fondue");
	doubleVisualDepList.set("Cheese to Fondue", "Swiss", "The cheese to fondue");
	double (*fondueFunc)(double);
	fondueFunc = fondueTempTester;

	Teuchos::RCP<Optika::NumberVisualDependency<double> > fondueDep = 
	Teuchos::RCP<Optika::NumberVisualDependency<double> >(
		new Optika::NumberVisualDependency<double>(
			"Temperature", 
			Teuchos::sublist(My_deplist,"NumberVisual Dependency List (double)"),
			"Cheese to Fondue", 
			Teuchos::sublist(My_deplist,"NumberVisual Dependency List (double)"),
			fondueFunc
		)
	);
	depSheet1->addDependency(fondueDep);
	fondueDep->evaluate();
	TEST_ASSERT(fondueDep->isDependentVisible());
	doubleVisualDepList.set("Temperature",99.0);
	fondueDep->evaluate();
	TEST_ASSERT(!fondueDep->isDependentVisible());

	/*
	 * Testing the BoolVisualDependency
	 */
	Teuchos::ParameterList&
	boolVisDepList = My_deplist->sublist(
		"Bool Visual Dependency List", 
		false,
		"Bool Visual Dependency testing list.\nWorking June 29 2009"
	);
	boolVisDepList.set("ShowPrecs", true, "Whether or not to should the Preciondtioner list");
	Teuchos::ParameterList&
	Prec_List0 = boolVisDepList.sublist("Preconditioner",false,"Sublist that defines the preconditioner.");
	Prec_List0.set("Type", "ILU", "The tpye of preconditioner to use");
	Teuchos::RCP<Optika::EnhancedNumberValidator<double> > droptolValidator = Teuchos::rcp(new Optika::EnhancedNumberValidator<double>(0,10,1e-3));
	Prec_List0.set("Drop Tolerance", 1e-3,"The tolerance below which entries from the\n""factorization are left out of the factors.", droptolValidator);
	Teuchos::RCP<Optika::BoolVisualDependency> 
	precDep1 = Teuchos::RCP<Optika::BoolVisualDependency>(
		new Optika::BoolVisualDependency(
			"ShowPrecs", 
			Teuchos::sublist(My_deplist,"Bool Visual Dependency List"),
			"Preconditioner", 
			Teuchos::sublist(My_deplist,"Bool Visual Dependency List"),
			true
		)
	);
	depSheet1->addDependency(precDep1);
	precDep1->evaluate();
	TEST_ASSERT(precDep1->isDependentVisible());
	boolVisDepList.set("ShowPrecs", false);
	precDep1->evaluate();
	TEST_ASSERT(!precDep1->isDependentVisible());



	/*
	 * Testing the StringVisualDepenency
	 */
	Teuchos::ParameterList&
    stringVisDepList = My_deplist->sublist(
		"String Visual Dependency List",
		false,
		"String Visual Dependency testing list.\nWorking June 29 2009"
	);
	Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
	favCheeseValidator = Teuchos::rcp(
		new Teuchos::StringToIntegralParameterEntryValidator<int>(
			Teuchos::tuple<std::string>( "Swiss", "American", "Cheder" ),
			"Favorite Cheese"
		)
	);
   
	stringVisDepList.set("Favorite Cheese", "American", "Your favorite type of cheese", favCheeseValidator);
	Teuchos::RCP<Optika::EnhancedNumberValidator<int> > 
	swissValidator = Teuchos::rcp(new Optika::EnhancedNumberValidator<int>(0,10));
	stringVisDepList.set("Swiss rating", 0, "How you rate swiss on a scale of 1 to 10", swissValidator);
	Teuchos::RCP<Optika::StringVisualDependency> 
	swissDep1 = Teuchos::RCP<Optika::StringVisualDependency>(
		new Optika::StringVisualDependency(
			"Favorite Cheese", 
			Teuchos::sublist(My_deplist,"String Visual Dependency List"),
			"Swiss rating", 
			Teuchos::sublist(My_deplist,"String Visual Dependency List"),
			"Swiss", 
			true
		)
	);
	depSheet1->addDependency(swissDep1);
	swissDep1->evaluate();
	TEST_ASSERT(!swissDep1->isDependentVisible());
	stringVisDepList.set("Favorite Cheese", "Swiss");
	swissDep1->evaluate();
	TEST_ASSERT(swissDep1->isDependentVisible());

	/*
	 * String Visual Tester with multiple values
	 */
	Teuchos::ParameterList&
    multiStringVisDepList = My_deplist->sublist(
		"Multi String Visual Dependency List",
		false
	);
	Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
	favCheeseValidator2 = Teuchos::rcp(
		new Teuchos::StringToIntegralParameterEntryValidator<int>(
			Teuchos::tuple<std::string>( "Provalone", "Swiss", "American", "Cheder" ),
			"Favorite Cheese"
		)
	);
   
	multiStringVisDepList.set("Favorite Cheese", "American", "Your favorite type of cheese", favCheeseValidator2);
	multiStringVisDepList.set("Swiss rating", 0, "How you rate swiss on a scale of 1 to 10", swissValidator);
	Teuchos::RCP<Optika::StringVisualDependency> 
	swissDep2 = Teuchos::RCP<Optika::StringVisualDependency>(
		new Optika::StringVisualDependency(
			"Favorite Cheese", 
			Teuchos::sublist(My_deplist,"Multi String Visual Dependency List"),
			"Swiss rating", 
			Teuchos::sublist(My_deplist,"Multi String Visual Dependency List"),
			Teuchos::tuple<std::string>("Swiss", "Cheder"), 
			true
		)
	);
	depSheet1->addDependency(swissDep2);
	swissDep2->evaluate();
	TEST_ASSERT(!swissDep2->isDependentVisible());
	multiStringVisDepList.set("Favorite Cheese", "Cheder");
	swissDep2->evaluate();
	TEST_ASSERT(swissDep2->isDependentVisible());

	/*
	 * Another test of the NumberVisualDependency.
	 */
	int (*visfunc)(int);
	visfunc = intVisualTester;
	Teuchos::ParameterList&
    numberVisDepList = My_deplist->sublist(
		"Number Visual Dependency List", 
		false, 
		"Number Visual Dependency testing list.\nWorking June 27th 2009"
	);
	numberVisDepList.set("Ice", 50, "Ice stuff");
	numberVisDepList.set("Room Temp", 10, "Room temperature");
	Teuchos::RCP<Optika::NumberVisualDependency<int> > 
	iceDep = Teuchos::RCP<Optika::NumberVisualDependency<int> >(
		new Optika::NumberVisualDependency<int>(
			"Room Temp", 
			Teuchos::sublist(My_deplist,"Number Visual Dependency List"),
			"Ice", 
			Teuchos::sublist(My_deplist,"Number Visual Dependency List"),
			visfunc
		)
	);
	depSheet1->addDependency(iceDep);
	iceDep->evaluate();
	TEST_ASSERT(iceDep->isDependentVisible());
	numberVisDepList.set("Room Temp", 33);
	iceDep->evaluate();
	TEST_ASSERT(!iceDep->isDependentVisible());

	/*
	 * Test condition visual dependency
	 */
	Teuchos::RCP<Teuchos::ParameterList> conVisDepList = Teuchos::sublist(My_deplist,"Condition Visual Dependency List", false);
	conVisDepList->set("double param", 4.0, "double parameter");
	conVisDepList->set("bool param", true, "bool parameter");
	conVisDepList->set("string param", "blah", "a string parameter");
	Teuchos::RCP<Optika::NumberCondition<double> > numberCon = Teuchos::rcp( new Optika::NumberCondition<double>("double param", conVisDepList, true));
	Teuchos::RCP<Optika::BoolCondition> boolCon = Teuchos::rcp(new Optika::BoolCondition("bool param", conVisDepList));
	Optika::Condition::ConditionList conList = Teuchos::tuple<Teuchos::RCP<Optika::Condition> >(numberCon, boolCon);
	Teuchos::RCP<Optika::AndCondition> andCon = Teuchos::rcp(new Optika::AndCondition(conList));
	Teuchos::RCP<Optika::ConditionVisualDependency> conVisDep = Teuchos::rcp(new Optika::ConditionVisualDependency(andCon, "string param", conVisDepList, true));
	depSheet1->addDependency(conVisDep);
	conVisDep->evaluate();
	TEST_ASSERT(conVisDep->isDependentVisible());
	conVisDepList->set("bool param", false);
	conVisDep->evaluate();
	TEST_ASSERT(!conVisDep->isDependentVisible());



	return (success ? 0:1);
}


/**
 * Test the ArrayLengthDependency.
 */
int testArrayLengthDep(Teuchos::FancyOStream &out){
	bool success = true;
	Teuchos::RCP<Teuchos::ParameterList> My_deplist = Teuchos::RCP<Teuchos::ParameterList>(new Teuchos::ParameterList);
	Teuchos::RCP<Optika::DependencySheet> depSheet1 = Teuchos::RCP<Optika::DependencySheet>(new Optika::DependencySheet(My_deplist));

	Teuchos::ParameterList&
	numberArrayLengthDepList = My_deplist->sublist("Number Array Length Dependency List", false, "Number Array Length Dependecy testing list.\nWorking June 27th 2009");
	numberArrayLengthDepList.set("Array Length", 10, "array length setter");
	Teuchos::Array<double> variableLengthArray(10,23.0);
	Teuchos::RCP<Optika::EnhancedNumberValidator<double> > 
	varLengthArrayVali = Teuchos::RCP<Optika::EnhancedNumberValidator<double> >(
  		new Optika::EnhancedNumberValidator<double>(10,50,4) 
	);
	numberArrayLengthDepList.set("Variable Length Array", variableLengthArray, "variable length array",
	Teuchos::RCP<Optika::ArrayNumberValidator<double> >(new Optika::ArrayNumberValidator<double>(varLengthArrayVali)));

	Teuchos::RCP<Optika::NumberArrayLengthDependency> 
	arrayLengthDep = Teuchos::RCP<Optika::NumberArrayLengthDependency>(
  		new Optika::NumberArrayLengthDependency(
			"Array Length", 
			Teuchos::sublist(My_deplist,"Number Array Length Dependency List"),
			"Variable Length Array", 
			Teuchos::sublist(My_deplist,"Number Array Length Dependency List")
		)
	);
	depSheet1->addDependency(arrayLengthDep);
	Teuchos::Array<double> dummyType;
	TEST_ASSERT(numberArrayLengthDepList.get("Variable Length Array", dummyType).length() ==10);
	numberArrayLengthDepList.set("Array Length", 12);
	arrayLengthDep()->evaluate();
	TEST_ASSERT(numberArrayLengthDepList.get("Variable Length Array", dummyType).length() ==12);
	numberArrayLengthDepList.set("Array Length", -1);
	TEST_THROW(arrayLengthDep()->evaluate(), Teuchos::Exceptions::InvalidParameterValue);

	return (success ? 0:1);
}

/**
 * Test the NumberValidatorAspectDependency.
 */
int testNumberValiAspDep(Teuchos::FancyOStream &out){
	bool success = true;
	Teuchos::RCP<Teuchos::ParameterList> My_deplist = Teuchos::RCP<Teuchos::ParameterList>(new Teuchos::ParameterList);
	Teuchos::RCP<Optika::DependencySheet> depSheet1 = Teuchos::RCP<Optika::DependencySheet>(new Optika::DependencySheet(My_deplist));

	Teuchos::ParameterList&
	numberValiAspDepList = My_deplist->sublist(
		"Number Validator Aspect Dependency List",
		false,
		"Number Validator Aspect Dependency testing list.\nWorking June 27th 2009"
	);
	Teuchos::RCP<Optika::EnhancedNumberValidator<int> > intVali2 = 
	Teuchos::rcp(new Optika::EnhancedNumberValidator<int>(0,20));
	numberValiAspDepList.set("Int", 8, "Int tester", intVali2);
	numberValiAspDepList.set("Int2", 8, "int2 tester", intVali2);
	numberValiAspDepList.set("Int dependee", 1, "Int dependee");

	int (*func)(int);
	func = intFuncTester;

	Teuchos::RCP<Optika::NumberValidatorAspectDependency<int> > 
	intDep1 = Teuchos::RCP<Optika::NumberValidatorAspectDependency<int> >(
		new Optika::NumberValidatorAspectDependency<int>(
			"Int dependee",
			Teuchos::sublist(My_deplist,"Number Validator Aspect Dependency List"),
			"Int",
			Teuchos::sublist(My_deplist,"Number Validator Aspect Dependency List"),
			intVali2,
			Optika::NumberValidatorAspectDependency<int>::Max,
			func
		)
	);
	
	Teuchos::RCP<Optika::NumberValidatorAspectDependency<int> > 
	intDep2 = Teuchos::RCP<Optika::NumberValidatorAspectDependency<int> >(
		new Optika::NumberValidatorAspectDependency<int>(
			"Int dependee",
			Teuchos::sublist(My_deplist,"Number Validator Aspect Dependency List"),
			"Int2",
			Teuchos::sublist(My_deplist,"Number Validator Aspect Dependency List"),
			intVali2,
			Optika::NumberValidatorAspectDependency<int>::Max,
			func
		)
	);

	depSheet1->addDependency(intDep1);
	depSheet1->addDependency(intDep2);

	TEST_ASSERT(Teuchos::rcp_static_cast<const Optika::EnhancedNumberValidator<int> >(numberValiAspDepList.getEntry("Int").validator())->max() == 20);
	TEST_ASSERT(Teuchos::rcp_static_cast<const Optika::EnhancedNumberValidator<int> >(numberValiAspDepList.getEntry("Int2").validator())->max() == 20);
	intDep1->evaluate();
	TEST_ASSERT(Teuchos::rcp_static_cast<const Optika::EnhancedNumberValidator<int> >(numberValiAspDepList.getEntry("Int").validator())->max() == 11);
	TEST_ASSERT(Teuchos::rcp_static_cast<const Optika::EnhancedNumberValidator<int> >(numberValiAspDepList.getEntry("Int2").validator())->max() == 11);


	return (success ? 0:1);
}

/**
 * Tests the excpetions associated with Dependencies
 */
int testDepExceptions(Teuchos::FancyOStream &out){
	bool success = true;
	Teuchos::RCP<Teuchos::ParameterList> list1 = Teuchos::RCP<Teuchos::ParameterList>(new Teuchos::ParameterList());
	Teuchos::RCP<Teuchos::ParameterList> list2 = Teuchos::RCP<Teuchos::ParameterList>(new Teuchos::ParameterList());

	list1->set("int parameter", 4, "int parameter");
	list1->set("double parameter", 6.0, "double parameter");
	list1->set("string parameter", "hahahaha", "string parameter");
	Teuchos::Array<double> doubleArray(10,23.0);
	list1->set("array parameter", doubleArray, "array parameter");
	list1->set("bool parameter", true, "bool parameter");

	/*
	 * Testing StringVisualDepenendcy exceptions.
	 */
	Teuchos::RCP<Optika::StringVisualDependency> stringVisDep;
	TEST_THROW(stringVisDep = Teuchos::RCP<Optika::StringVisualDependency>(new Optika::StringVisualDependency("not in list", list1, "int parameter", list1, "cheese", true)), Optika::InvalidDependencyException);
	TEST_THROW(stringVisDep = Teuchos::RCP<Optika::StringVisualDependency>(new Optika::StringVisualDependency("string parameter", list1, "not in list", list1, "cheese", true)), Optika::InvalidDependencyException);
	TEST_THROW(stringVisDep = Teuchos::RCP<Optika::StringVisualDependency>(new Optika::StringVisualDependency("double parameter", list1, "int parameter", list1, "cheese", true)), Optika::InvalidDependencyException);

	TEST_THROW(Teuchos::RCP<Optika::BoolVisualDependency> boolVisDep = Teuchos::RCP<Optika::BoolVisualDependency>(new Optika::BoolVisualDependency("int parameter", list1, "double parameter", list1, false)), Optika::InvalidDependencyException);

	TEST_THROW(Teuchos::RCP<Optika::NumberArrayLengthDependency> numArrayLengthDep = Teuchos::RCP<Optika::NumberArrayLengthDependency>(new Optika::NumberArrayLengthDependency("double parameter", list1, "array parameter", list1)), Optika::InvalidDependencyException);
	TEST_THROW(Teuchos::RCP<Optika::NumberArrayLengthDependency> numArrayLengthDep = Teuchos::RCP<Optika::NumberArrayLengthDependency>(new Optika::NumberArrayLengthDependency("int parameter", list1, "double parameter", list1)), Optika::InvalidDependencyException);

	/*
	 * Testing StringValidatorDependency exceptions.
	 */
	Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
    cheeseValidator = Teuchos::rcp(
		new Teuchos::StringToIntegralParameterEntryValidator<int>(
   			Teuchos::tuple<std::string>( "Swiss", "American", "Super Awesome Cheese" )
			,"Food Selector"
			)
	);

	Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
	sodaValidator = Teuchos::rcp(
		new Teuchos::StringToIntegralParameterEntryValidator<int>(
			Teuchos::tuple<std::string>( "Pepsi", "Coke", "Kurtis Cola", "Bad Cola" )
			,"Food Selector"
			)
		);

	Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
	chipsValidator = Teuchos::rcp(
		new Teuchos::StringToIntegralParameterEntryValidator<int>(
			Teuchos::tuple<std::string>( "Lays", "Doritos", "Kurtis Super Awesome Brand" )
			,"Food Selector"
		)
	);


	list1->set("string 2 parameter", "Swiss", "second string parameter", cheeseValidator);
	Optika::StringValidatorDependency::ValueToValidatorMap testValidatorMap1;
	testValidatorMap1["Cheese"] = cheeseValidator;
	testValidatorMap1["Soda"] = sodaValidator;
	testValidatorMap1["Chips"] = chipsValidator;
	TEST_THROW(Teuchos::RCP<Optika::StringValidatorDependency> stringValiDep = Teuchos::RCP<Optika::StringValidatorDependency>(new Optika::StringValidatorDependency("int parameter", list1, "string 2 parameter", list1, testValidatorMap1, cheeseValidator)), Optika::InvalidDependencyException);
	Teuchos::RCP<Optika::EnhancedNumberValidator<int> > intVali = Teuchos::rcp(new Optika::EnhancedNumberValidator<int>(0,20));
	testValidatorMap1["Candy"] = intVali;
	TEST_THROW(Teuchos::RCP<Optika::StringValidatorDependency> stringValiDep = Teuchos::RCP<Optika::StringValidatorDependency>(new Optika::StringValidatorDependency("string parameter", list1, "string 2 parameter", list1, testValidatorMap1, cheeseValidator)), Optika::InvalidDependencyException);
	
	/*
	 * Testing BoolValidatorDependency exceptions.
	 */
	Teuchos::RCP<Optika::EnhancedNumberValidator<double> > doubleVali1 = Teuchos::rcp(new Optika::EnhancedNumberValidator<double>(0.0,20.0));
	Teuchos::RCP<Optika::EnhancedNumberValidator<double> > doubleVali2 = Teuchos::rcp(new Optika::EnhancedNumberValidator<double>(5.0,20.0));
	list1->set("double parameter", 6.0, "double parameter", doubleVali1);
	TEST_THROW(Teuchos::RCP<Optika::BoolValidatorDependency> boolValiDep = Teuchos::RCP<Optika::BoolValidatorDependency>(new Optika::BoolValidatorDependency("int parameter", list1, "double parameter", list1, doubleVali1, doubleVali2)), Optika::InvalidDependencyException);
	TEST_THROW(Teuchos::RCP<Optika::BoolValidatorDependency> boolValiDep = Teuchos::RCP<Optika::BoolValidatorDependency>(new Optika::BoolValidatorDependency("bool parameter", list1, "double parameter", list1, intVali, doubleVali2)), Optika::InvalidDependencyException);
	TEST_THROW(Teuchos::RCP<Optika::BoolValidatorDependency> boolValiDep = Teuchos::RCP<Optika::BoolValidatorDependency>(new Optika::BoolValidatorDependency("bool parameter", list1, "double parameter", list1, doubleVali1, intVali)), Optika::InvalidDependencyException);

	TEST_THROW(Teuchos::RCP<Optika::NumberVisualDependency<int> > boolValiDep = Teuchos::RCP<Optika::NumberVisualDependency<int> >(new Optika::NumberVisualDependency<int>("bool parameter", list1, "double parameter", list1, intFuncTester)), Optika::InvalidDependencyException);


	/*
	 * Testing NumberValidatorAspectDependency exceptions.
	 */
	list1->set("int 2 parameter", 6, "int 2 parameter");
	TEST_THROW(Teuchos::RCP<Optika::NumberValidatorAspectDependency<int> > 
		intDep1 = Teuchos::RCP<Optika::NumberValidatorAspectDependency<int> >(
			new Optika::NumberValidatorAspectDependency<int>(
				"int parameter",
				list1,
				"int 2 parameter",
				list1,
				intVali,
				Optika::NumberValidatorAspectDependency<int>::Max,
				intFuncTester
			)
		),
		Optika::InvalidDependencyException
	);

	list1->set("int 2 parameter", 6, "int 2 parameter", intVali);
	TEST_THROW(Teuchos::RCP<Optika::NumberValidatorAspectDependency<int> > 
		intDep1 = Teuchos::RCP<Optika::NumberValidatorAspectDependency<int> >(
			new Optika::NumberValidatorAspectDependency<int>(
				"string parameter",
				list1,
				"int 2 parameter",
				list1,
				intVali,
				Optika::NumberValidatorAspectDependency<int>::Max,
				intFuncTester
			)
		),
		Optika::InvalidDependencyException
	);
	Teuchos::RCP<Optika::EnhancedNumberValidator<int> > intVali2 = Teuchos::rcp(new Optika::EnhancedNumberValidator<int>(5,20));
	TEST_THROW(Teuchos::RCP<Optika::NumberValidatorAspectDependency<int> > 
		intDep1 = Teuchos::RCP<Optika::NumberValidatorAspectDependency<int> >(
			new Optika::NumberValidatorAspectDependency<int>(
				"int parameter",
				list1,
				"int 2 parameter",
				list1,
				intVali2,
				Optika::NumberValidatorAspectDependency<int>::Max,
				intFuncTester
			)
		),
		Optika::InvalidDependencyException
	);

	TEST_THROW(Teuchos::RCP<Optika::NumberValidatorAspectDependency<int> > 
		intDep1 = Teuchos::RCP<Optika::NumberValidatorAspectDependency<int> >(
			new Optika::NumberValidatorAspectDependency<int>(
				"int parameter",
				list1,
				"double parameter",
				list1,
				intVali,
				Optika::NumberValidatorAspectDependency<int>::Max,
				intFuncTester
			)
		),
		Optika::InvalidDependencyException
	);

	TEST_THROW(Teuchos::RCP<Optika::NumberValidatorAspectDependency<int> > 
		intDep1 = Teuchos::RCP<Optika::NumberValidatorAspectDependency<int> >(
			new Optika::NumberValidatorAspectDependency<int>(
				"double parameter",
				list1,
				"int parameter",
				list1,
				intVali,
				Optika::NumberValidatorAspectDependency<int>::Max,
				intFuncTester
			)
		),
		Optika::InvalidDependencyException
	);

	/*
	 * Testing RangeValidatorDependency exceptions.
	 */
	list1->set("Cheese to Fondue", "Swiss", "the cheese to fondue");
	Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
	lowTempCheeseValidator = Teuchos::rcp(
		new Teuchos::StringToIntegralParameterEntryValidator<int>(
			Teuchos::tuple<std::string>( "PepperJack", "Swiss", "American" ),
			"Cheese to Fondue"
		)
	);
	Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
	highTempCheeseValidator = Teuchos::rcp(
		new Teuchos::StringToIntegralParameterEntryValidator<int>(
			Teuchos::tuple<std::string>( "Munster", "Provalone", "Kurtis Super Awesome Cheese"),
			"Cheese to Fondue"
		)
	);
	list1->set("Cheese to Fondue", "Swiss", "the cheese to fondue", lowTempCheeseValidator);
	Optika::RangeValidatorDependency<double>::RangeToValidatorMap tempranges;
	tempranges[std::pair<double,double>(100,200)] = lowTempCheeseValidator;
	tempranges[std::pair<double,double>(200,300)] = highTempCheeseValidator;
	TEST_THROW(
		Teuchos::RCP<Optika::RangeValidatorDependency<double> > 
		cheeseTempDep = Teuchos::RCP<Optika::RangeValidatorDependency<double> >(
			new Optika::RangeValidatorDependency<double>(
				"string parameter", 
				list1,	
				"Cheese to Fondue", 
				list1,	
				tempranges, 
				lowTempCheeseValidator
			)
		),
		Optika::InvalidDependencyException
	);
	tempranges[std::pair<double,double>(400,800)] = intVali;
	TEST_THROW(
		Teuchos::RCP<Optika::RangeValidatorDependency<double> > 
		cheeseTempDep = Teuchos::RCP<Optika::RangeValidatorDependency<double> >(
			new Optika::RangeValidatorDependency<double>(
				"int parameter", 
				list1,	
				"Cheese to Fondue", 
				list1,
				tempranges, 
				lowTempCheeseValidator
			)
		),
		Optika::InvalidDependencyException
	);


	/*
	 * Testing DependencySheet exceptions.
	 */
	Teuchos::RCP<Optika::DependencySheet> depSheet1 = Teuchos::RCP<Optika::DependencySheet>(new Optika::DependencySheet(list1));
	Teuchos::RCP<Optika::DependencySheet> depSheet2 = Teuchos::RCP<Optika::DependencySheet>(new Optika::DependencySheet(list2));

	list2->set("list2 double", 4.0, "a double parameter in list 2");
	list2->set("list2 bool", true, "a bool parameter in list2");

	Teuchos::RCP<Optika::BoolVisualDependency> boolVisDep = Teuchos::RCP<Optika::BoolVisualDependency>(new Optika::BoolVisualDependency("bool parameter", list1, "list2 double", list2, false));
	TEST_THROW(depSheet2->addDependency(boolVisDep), Optika::InvalidDependencyException);
	boolVisDep = Teuchos::RCP<Optika::BoolVisualDependency>(new Optika::BoolVisualDependency("list2 bool", list2, "double parameter", list1, false));
	TEST_THROW(depSheet2->addDependency(boolVisDep), Optika::InvalidDependencyException);


	return (success ? 0:1);
}

int main(int argc, char* argv[]){
	bool success = true;
	Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
	if(testValiDeps(*out) == 1){
		success = false;
	}
	if(testVisualDeps(*out) == 1){
		success = false;
	}
	if(testArrayLengthDep(*out) == 1){
		success = false;
	}
	if(testNumberValiAspDep(*out) == 1){
		success = false;
	}

	if(testDepExceptions(*out) ==1){
		success = false;
	}
	return (success ? 0:1);
}

