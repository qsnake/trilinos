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
#ifndef OPTIKA_VALIDATORFACTORY_HPP_
#define OPTIKA_VALIDATORFACTORY_HPP_
#include "Optika_SpecificParameterEntryValidators.hpp"


namespace Optika{

class ValidatorFactory{
public:
	enum ValidatorType{Int, Short, Double, Float, FileName,
	IntArray, ShortArray, DoubleArray, FloatArray, FileNameArray};

	/**
	 * Creates a validator of the given type.
	 * 
	 * @param valiType The type of validator to be created.
	 * @return A validator of the specified type.
	 */
	static Teuchos::RCP<Teuchos::ParameterEntryValidator> createValidator(ValidatorType valiType);

	/**
	 * Creates and returns a Enhanced Number Validator of type int.
	 *
	 * @return An Enhanced Number Validator of type int.
	 */
	static Teuchos::RCP<EnhancedNumberValidator<int> > getIntValidator();

	/**
	 * Creates and returns a Enhanced Number Validator of type short.
	 *
	 * @return An Enhanced Number Validator of type short.
	 */
	static Teuchos::RCP<EnhancedNumberValidator<short> > getShortValidator();

	/**
	 * Creates and returns a Enhanced Number Validator of type double.
	 *
	 * @return An Enhanced Number Validator of type double.
	 */
	static Teuchos::RCP<EnhancedNumberValidator<double> > getDoubleValidator();

	/**
	 * Creates and returns a Enhanced Number Validator of type float.
	 *
	 * @return An Enhanced Number Validator of type float.
	 */
	static Teuchos::RCP<EnhancedNumberValidator<float> > getFloatValidator();

	/**
	 * Creates and returns FileNameValidator.
	 *
	 * @return A FileNameValidator.
	 */
	static Teuchos::RCP<FileNameValidator> getFileNameValidator();

	/**
	 * Creates and returns an Array Number Validator of type int.
	 *
	 * @return An Enhanced Number Validator of type int.
	 */
	static Teuchos::RCP<ArrayNumberValidator<int> > getArrayIntValidator();

	/**
	 * Creates and returns an Array Number Validator of type short.
	 *
	 * @return An Enhanced Number Validator of type short.
	 */
	static Teuchos::RCP<ArrayNumberValidator<short> > getArrayShortValidator();

	/**
	 * Creates and returns an Array Number Validator of type double.
	 *
	 * @return An Enhanced Number Validator of type double.
	 */
	static Teuchos::RCP<ArrayNumberValidator<double> > getArrayDoubleValidator();

	/**
	 * Creates and returns an Array Number Validator of type float.
	 *
	 * @return An Enhanced Number Validator of type float.
	 */
	static Teuchos::RCP<ArrayNumberValidator<float> > getArrayFloatValidator();

	/**
	 * Creates and returns an Array File Name Validator.
	 *
	 * @return An Array File Name Validator.
	 */
	static Teuchos::RCP<ArrayFileNameValidator> getArrayFileNameValidator();
};

}

#endif /* OPTIKA_VALIDATORFACTORY_HPP_ */
