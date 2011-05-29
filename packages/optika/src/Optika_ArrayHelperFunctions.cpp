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
#include "Optika_ArrayHelperFunctions.hpp"

namespace Optika{
bool doesParameterContainArray(const Teuchos::ParameterEntry *parameter){
	QString typeName = QString::fromStdString(parameter->getAny(false).typeName());
	return typeName.contains("Teuchos") && typeName.contains("Array");	
}

QStringList getValues(QString& values){
	values = values.remove("{");
	values = values.remove("}");
	QStringList toReturn = values.split(",");
	for(int i = 0; i < toReturn.size(); ++i){
		if(toReturn[i].at(0) == QChar(' ')){
			toReturn[i] = toReturn[i].remove(0,1);
		}
	}
	return toReturn;
}

QString determineArrayType(Teuchos::ParameterEntry *parameter){
	Teuchos::any anyArray = parameter->getAny();
	if(anyArray.type() == typeid(Teuchos::Array<int>)){
		return intId;
	}
	else if(anyArray.type() == typeid(Teuchos::Array<short>)){
		return shortId;
	}
	else if(anyArray.type() == typeid(Teuchos::Array<double>)){
		return doubleId;
	}
	else if(anyArray.type() == typeid(Teuchos::Array<float>)){
		return floatId;
	}
	else if(anyArray.type() == typeid(Teuchos::Array<std::string>)){
		return stringId;
	}
	else{
		return unrecognizedId;		
	}
}

template <>
Teuchos::Array<std::string> fromStringToArray<std::string>(QString arrayString){
	arrayString = arrayString.remove("{");
	arrayString = arrayString.remove("}");
	QStringList tempValues = arrayString.split(",");
	for(int i = 0; i < tempValues.size(); ++i){
		if(tempValues[i].at(0) == QChar(' ')){
			tempValues[i] = tempValues[i].remove(0,1);
		}
	}
	QList<QVariant> values;
	for(int i = 0; i<tempValues.size(); ++i){
		values.append(tempValues[i]);
	}
	Teuchos::Array<std::string> toReturn;
	for(int i = 0; i<values.size(); ++i){
		toReturn.append(values[i].value<QString>().toStdString());	
	}
	return toReturn;

}


}

