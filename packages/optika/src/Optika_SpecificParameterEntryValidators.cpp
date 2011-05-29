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
#include "Optika_SpecificParameterEntryValidators.hpp"

namespace Optika{

FileNameValidator::FileNameValidator(bool mustAlreadyExist):ParameterEntryValidator(),mustAlreadyExist(mustAlreadyExist){}

bool FileNameValidator::fileMustExist() const{
	return mustAlreadyExist;
}

bool FileNameValidator::setFileMustExist(bool shouldFileExist){
	this->mustAlreadyExist = shouldFileExist;
	return mustAlreadyExist;
}

Teuchos::RCP<const Teuchos::Array<std::string> > FileNameValidator::validStringValues() const{
	return Teuchos::null;
}

void FileNameValidator::validate(Teuchos::ParameterEntry const &entry, std::string const &paramName, std::string const &sublistName) const{
	Teuchos::any anyValue = entry.getAny(true);
	if(!(anyValue.type() == typeid(std::string) )){
		const std::string &entryName = entry.getAny(false).typeName();
		std::stringstream oss;
		std::string msg;
		oss << "Aww shoot! Sorry bud, but it looks like the \"" << paramName << "\"" <<
		" parameter in the \"" << sublistName << "\" sublist didn't quite work out.\n" <<
		"No need to fret though. I'm sure it's just a small mistake. Maybe the information below "<<
		"can help you figure out what went wrong.\n\n"
		"Error: The value that you entered was the wrong type.\n" <<
		"Parameter: " << paramName << "\n" << 
		"Type specified: " << entryName << "\n" <<
		"Type accepted: " << typeid(std::string).name() << "\n";
		msg = oss.str();
		throw Teuchos::Exceptions::InvalidParameterType(msg);
	}
	if(mustAlreadyExist){
		std::string fileName = Teuchos::getValue<std::string>(entry);
		struct stat fileInfo;
		int intStat= stat(fileName.c_str(),&fileInfo);
		if(intStat !=0){
			std::stringstream oss;
			std::string msg;
			oss << "Aww shoot! Sorry bud, but it looks like the \"" << paramName << "\"" <<
			" parameter in the \"" << sublistName << "\" sublist didn't quite work out.\n" <<
			"No need to fret though. I'm sure it's just a small mistake. Maybe the information below "<<
			"can help you figure out what went wrong.\n\n"
			"Error: The file must already exists. The value you entered does not corresspond to an existing file name.\n" <<
			"Parameter: " << paramName << "\n" << 
			"File name specified: " << fileName << "\n";
			msg = oss.str();
			throw Teuchos::Exceptions::InvalidParameterValue(msg);
		}
	}
}

void FileNameValidator::printDoc(std::string const &docString, std::ostream &out) const{
	Teuchos::StrUtils::printLines(out,"# ",docString);
	out << "#  Validator Used: \n";
	out << "#	FileName Validator\n";
}

StringValidator::StringValidator(ValueList validStrings):
	ParameterEntryValidator(),
	validStrings(validStrings){}

const Teuchos::Array<std::string> StringValidator::setValidStrings(ValueList validStrings){
	this->validStrings = validStrings;
	return this->validStrings;
}

Teuchos::RCP<const Teuchos::Array<std::string> > StringValidator::validStringValues() const{
	Teuchos::RCP<const Teuchos::Array<std::string> > toReturn = Teuchos::RCP<const Teuchos::Array<std::string> >(new Teuchos::Array<std::string>(validStrings));
	return toReturn;
}

void StringValidator::validate(Teuchos::ParameterEntry const &entry, std::string const &paramName, std::string const &sublistName) const{
	Teuchos::any anyValue = entry.getAny(true);
	if(!(anyValue.type() == typeid(std::string) )){
		const std::string &entryName = entry.getAny(false).typeName();
		std::stringstream oss;
		std::string msg;
		oss << "Aww shoot! Sorry bud, but it looks like the \"" << paramName << "\"" <<
		" parameter in the \"" << sublistName << "\" sublist didn't quite work out.\n" <<
		"No need to fret though. I'm sure it's just a small mistake. Maybe the information below "<<
		"can help you figure out what went wrong.\n\n"
		"Error: The value that you entered was the wrong type." <<
		"Parameter: " << paramName << "\n" << 
		"Type specified: " << entryName << "\n" <<
		"Type accepted: " << typeid(std::string).name() << "\n";
		msg = oss.str();
		throw Teuchos::Exceptions::InvalidParameterType(msg);
	}
	else{
		Teuchos::Array<std::string>::const_iterator it = std::find(validStrings.begin(), validStrings.end(), Teuchos::getValue<std::string>(entry));
		if(it == validStrings.end()){
			std::stringstream oss;
			std::string msg;
			oss << "Aww shoot! Sorry bud, but it looks like the \"" << paramName << "\"" <<
			" parameter in the \"" << sublistName << "\" sublist didn't quite work out.\n" <<
			"No need to fret though. I'm sure it's just a small mistake. Maybe the information below "<<
			"can help you figure out what went wrong.\n\n"
			"Error: The value that was entered doesn't fall with in " <<
			"the range set by the validator." <<
			"Parameter: " << paramName << "\n" <<
			"Acceptable Values: " << validStrings << "\n" <<
			"Value entered: " << Teuchos::getValue<std::string>(entry) << "\n";
			msg = oss.str();
			throw Teuchos::Exceptions::InvalidParameterValue(msg);
		}

	}
}

void StringValidator::printDoc(std::string const &docString, std::ostream &out) const{
	Teuchos::StrUtils::printLines(out,"# ",docString);
	out << "#  Validator Used: \n";
	out << "#	FileName Validator\n";
}

}

