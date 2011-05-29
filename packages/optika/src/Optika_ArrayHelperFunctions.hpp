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
#ifndef OPTIKA_ARRAYHELPERFUNCTIONS_HPP_
#define OPTIKA_ARRAYHELPERFUNCTIONS_HPP_
#include <QStringList>
#include <QVariant>
#include "Optika_Types.hpp"
#include "Teuchos_ParameterEntry.hpp"
namespace Optika{

/**
 * Determines whether or not a ParameterEntry contains an array.
 *
 * @return True if the ParameterEntry contains an array, false otherwise.
 */
bool doesParameterContainArray(const Teuchos::ParameterEntry *parameter);

/**
 * Takes a string representing an array, formats it, and returns
 * a QStringList containing each value in the array.
 *
 * @param values A QString containing the values in the array.
 * @return A QStringList containing the values in the array.
 */
QStringList getValues(QString& values);

/**
 * Determines the type of array stored in a parameter.
 *
 * @param parameter The parameter whose array type is in question.
 * @return A QString containing the type of array in the parameter.
 */
QString determineArrayType(Teuchos::ParameterEntry *parameter);

template <class S>
Teuchos::Array<S> fromStringToArray(QString arrayString){
	arrayString = arrayString.remove("{");
	arrayString = arrayString.remove("}");
	arrayString = arrayString.remove(" ");
	QStringList tempValues = arrayString.split(",");
	QList<QVariant> values;
	for(int i = 0; i<tempValues.size(); ++i){
		values.append(tempValues[i]);
	}
	Teuchos::Array<S> toReturn;
	for(int i = 0; i<values.size(); ++i){
		toReturn.append(values[i].value<S>());	
	}
	return toReturn;

}

template <>
Teuchos::Array<std::string> fromStringToArray<std::string>(QString arrayString);

}
#endif /* OPTIKA_ARRAYHELPERFUNCTIONS_HPP_ */
