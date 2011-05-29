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
#include "Optika_ValidatorFactory.hpp"
#include "Teuchos_ENull.hpp"
namespace Optika{

Teuchos::RCP<Teuchos::ParameterEntryValidator> ValidatorFactory::createValidator(ValidatorType valiType){
	switch(valiType){
		case Int:
			return Teuchos::RCP<EnhancedNumberValidator<int> >(new EnhancedNumberValidator<int>());
			break;
		case Short:
			return Teuchos::RCP<EnhancedNumberValidator<short> >(new EnhancedNumberValidator<short>());
			break;
		case Double:
			return Teuchos::RCP<EnhancedNumberValidator<double> >(new EnhancedNumberValidator<double>());
			break;
		case Float:
			return Teuchos::RCP<EnhancedNumberValidator<float> >(new EnhancedNumberValidator<float>());
			break;
		case IntArray:
			return Teuchos::RCP<ArrayNumberValidator<int> >( new ArrayNumberValidator<int>( Teuchos::RCP<EnhancedNumberValidator<int> >( new EnhancedNumberValidator<int>())));
			break;
		case ShortArray:
			return Teuchos::RCP<ArrayNumberValidator<short> >( new ArrayNumberValidator<short>( Teuchos::RCP<EnhancedNumberValidator<short> >( new EnhancedNumberValidator<short>())));
			break;
		case DoubleArray:
			return Teuchos::RCP<ArrayNumberValidator<double> >( new ArrayNumberValidator<double>( Teuchos::RCP<EnhancedNumberValidator<double> >( new EnhancedNumberValidator<double>())));
			break;
		case FloatArray:
			return Teuchos::RCP<ArrayNumberValidator<float> >( new ArrayNumberValidator<float>(Teuchos::RCP<EnhancedNumberValidator<float> >( new EnhancedNumberValidator<float>())));
			break;
		case FileName:
			return Teuchos::RCP<FileNameValidator>(new FileNameValidator());
			break;
		case FileNameArray:
			return Teuchos::RCP<ArrayFileNameValidator>(new ArrayFileNameValidator(Teuchos::RCP<FileNameValidator>(new FileNameValidator())));
			break;
		default:
			Teuchos::RCP<Teuchos::ParameterEntryValidator> toReturn;
			return toReturn;
			break;
	}
	Teuchos::RCP<Teuchos::ParameterEntryValidator> toReturn;
	return toReturn;
}

Teuchos::RCP<EnhancedNumberValidator<int> > ValidatorFactory::getIntValidator(){
	return Teuchos::RCP<EnhancedNumberValidator<int> >(new EnhancedNumberValidator<int>());
}

Teuchos::RCP<EnhancedNumberValidator<short> > ValidatorFactory::getShortValidator(){
	return Teuchos::RCP<EnhancedNumberValidator<short> >(new EnhancedNumberValidator<short>());
}

Teuchos::RCP<EnhancedNumberValidator<double> > ValidatorFactory::getDoubleValidator(){
	return Teuchos::RCP<EnhancedNumberValidator<double> >(new EnhancedNumberValidator<double>());
}

Teuchos::RCP<EnhancedNumberValidator<float> > ValidatorFactory::getFloatValidator(){
	return Teuchos::RCP<EnhancedNumberValidator<float> >(new EnhancedNumberValidator<float>());
}

Teuchos::RCP<FileNameValidator> ValidatorFactory::getFileNameValidator(){
	return Teuchos::RCP<FileNameValidator>(new FileNameValidator());
}

Teuchos::RCP<ArrayNumberValidator<int> > ValidatorFactory::getArrayIntValidator(){
	return Teuchos::RCP<ArrayNumberValidator<int> >( new ArrayNumberValidator<int>( Teuchos::RCP<EnhancedNumberValidator<int> >( new EnhancedNumberValidator<int>())));
}

Teuchos::RCP<ArrayNumberValidator<short> > ValidatorFactory::getArrayShortValidator(){
	return Teuchos::RCP<ArrayNumberValidator<short> >( new ArrayNumberValidator<short>( Teuchos::RCP<EnhancedNumberValidator<short> >( new EnhancedNumberValidator<short>())));
}

Teuchos::RCP<ArrayNumberValidator<double> > ValidatorFactory::getArrayDoubleValidator(){
	return Teuchos::RCP<ArrayNumberValidator<double> >( new ArrayNumberValidator<double>( Teuchos::RCP<EnhancedNumberValidator<double> >( new EnhancedNumberValidator<double>())));
}

Teuchos::RCP<ArrayNumberValidator<float> > ValidatorFactory::getArrayFloatValidator(){
	return Teuchos::RCP<ArrayNumberValidator<float> >( new ArrayNumberValidator<float>( Teuchos::RCP<EnhancedNumberValidator<float> >( new EnhancedNumberValidator<float>())));
}

Teuchos::RCP<ArrayFileNameValidator> ValidatorFactory::getArrayFileNameValidator(){
	return Teuchos::RCP<ArrayFileNameValidator>(new ArrayFileNameValidator(Teuchos::RCP<FileNameValidator>(new FileNameValidator())));
}

}
