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
#ifndef OPTIKA_SPECIFICPARAMETERENTRYVALIDATORS_HPP_
#define OPTIKA_SPECIFICPARAMETERENTRYVALIDATORS_HPP_
#include <iostream>
#include <sstream>
#include <string>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <float.h>
#include <limits>
#include <sys/stat.h>
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Optika_ArrayHelperFunctions.hpp"

namespace Optika{


static int intDefaultStep = 1;
static short shortDefaultStep = 1;
//static long long int longlongDefaultStep = 1;
static double doubleDefaultStep = 1;
static int doubleDefaultPrecision = 3;
static float floatDefaultStep = 1;
static int floatDefaultPrecision = 3;

/**
 * A Template base class for NumberValidators.
 * Note that while this is not an abstract base class,
 * you really shouldn't use it.
 * Just use the EnhancedNumberValidator class.
 */
template <class S>
class GenericNumberValidator : public Teuchos::ParameterEntryValidator{
public:
	/**
	 * Constructs a GenericNumberValidator.
	 *
	 * @param type The type of the validator.
	 * @param min The minimum acceptable value for this validator.
	 * @param max The maximum acceptable value for this validator.
	 * @param step The increments at which the value should be changed. This is mostly used for 
	 * the QSpinBox that is used in the Optika GUI. If you're not using the GUI, you may ignore this parameter.
	 */
	GenericNumberValidator(QString type, S min, S max, S step):Teuchos::ParameterEntryValidator(),
	type(type), minVal(min), maxVal(max), step(step), containsMin(true), containsMax(true){}

	/**
	 * Constructs a GenericNumberValidator without an explicit minimum or maximum.
	 *
	 * @param type The type of the validator.
	 * @param step The increments at which the value should be changed. This is mostly used for 
	 * the QSpinBox that is used in the Optika GUI. If you're not using the GUI, you may ignore this parameter.
	 */
	GenericNumberValidator(QString type, S step):
		Teuchos::ParameterEntryValidator(),
		type(type),
		step(step),
		containsMin(false),
		containsMax(false)
	{
		if(std::numeric_limits<S>::is_integer){
			this->minVal = std::numeric_limits<S>::min();
			this->maxVal = std::numeric_limits<S>::max();
		}
		else{
			this->minVal = -std::numeric_limits<S>::max();
			this->maxVal = std::numeric_limits<S>::max();
		}
	}
		
	/**
	 * Sets the minimum acceptable value for the validator.
	 * 
	 * @param min The desired minimum acceptable value for the validator.
	 */
	void setMin(S min){
		minVal = min;
		containsMin = true;
	}

	/**
	 * Sets the maximum acceptable value for the validator.
	 * 
	 * @param min The desired maximum acceptable value for the validator.
	 */
	void setMax(S max){
		maxVal = max;
		containsMax = true;
	}

	/**
	 * Gets the minimum acceptable value for the validator.
	 *
	 *@return The minimum acceptable value for the validator.
	 */
	S min() const{
		return minVal;
	}

	/**
	 * Gets the maximum acceptable value for the validator.
	 *
	 *@return The maximum acceptable value for the validator.
	 */
	S max() const{
		return maxVal;
	}

	/**
	 * Determines whether or not the validator has a minimum value.
	 *
	 * @return True if the validator has a minimum value, false otherwise.
	 */
	bool hasMin() const{
		return containsMin;
	}

	/**
	 * Determines whether or not the validator has a maximum value.
	 *
	 * @return True if the validator has a maximum value, false otherwise.
	 */ 
	bool hasMax() const{
		return containsMax;
	}

	/**
	 * Gets the step being used for the validator.
	 *
	 * @return The step being used for the validator.
	 */
	S getStep() const{
		return step;
	}

	/**
	 * Sets the step being used for the validator.
	 *
	 * @param The step to be used for the validator.
	 */
	void setStep(S step){
		this->step = step;
	}

	/**
	 * Gets the type of number being used.
	 *
	 * @return A string containting the name of the type.
	 */
	virtual const QString getType() const{
		return type;
	}

	Teuchos::RCP< const Teuchos::Array<std::string> > validStringValues() const{
		return Teuchos::null;
	}

	void validate(Teuchos::ParameterEntry const &entry, std::string const &paramName, std::string const &sublistName) const{
		Teuchos::any anyValue = entry.getAny(true);
		if(anyValue.type() == typeid(S) ){
			if(!(Teuchos::any_cast<S>(anyValue) >= minVal && Teuchos::any_cast<S>(anyValue) <= maxVal)){
				std::stringstream oss;
				std::string msg;
				oss << "Aww shoot! Sorry bud, but it looks like the \"" << paramName << "\"" <<
				" parameter in the \"" << sublistName << "\" sublist didn't quite work out.\n" <<
				"No need to fret though. I'm sure it's just a small mistake. Maybe the information below "<<
				"can help you figure out what went wrong.\n\n"
				"Error: The value that was entered doesn't fall with in " <<
				"the range set by the validator.\n" <<
				"Parameter: " << paramName << "\n" <<
				"Min: " << minVal << "\n" <<
				"Max: " << maxVal << "\n" <<
				"Value entered: " << (Teuchos::any_cast<S>(anyValue)) << "\n";
				msg = oss.str();
				throw Teuchos::Exceptions::InvalidParameterValue(msg);
			}	
		}
		else{
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
			"Type accepted: " << typeid(S).name() << "\n";
			msg = oss.str();
			throw Teuchos::Exceptions::InvalidParameterType(msg);
		}
	}

	void printDoc(std::string const &docString, std::ostream &out) const{
		Teuchos::StrUtils::printLines(out,"# ",docString);
		out << "#  Validator Used: \n";
		out << "#  	Number Validator\n";
		out << "#  	Type: " << getType().toStdString() << "\n";
		out << "#  	Min (inclusive): " << minVal << "\n";
		out << "#  	Max (inclusive): " << maxVal << "\n";
	}

private:
	/**
	 * The type of the array.
	 */
	QString type;

	/**
	 * The minimum value accepted by the validator.
	 */
	S minVal;

	/**
	 * The maximum value accepted by the validator.
	 */
	S maxVal;

	/**
	 * The increment to use when increaseing or decreaseing the value the validator is validating.
	 */
	S step;

	/**
	 * Whether or not a minimum value has been specified for this validator.
	 */
	bool containsMin;

	/**
	 * Whetehr or not a maximum value has been specified for this validator.
	 */
	bool containsMax;
};


/**
 * Validates inputs using mins and max for specific number types.
 */
template <class S> 
class EnhancedNumberValidator : public GenericNumberValidator<S>{
public:
	/**
	 * Constructs an EnhancedNumberValidator.
	 *
	 * @param min The minimum acceptable value for this validator.
	 * @param max The maximum acceptable value for this validator.
	 * @param step The increments at which the value should be changed. This is mostly used for 
	 * the QSpinBox that is used in the Optika GUI. If you're not using the GUI, you may ignore this parameter.
	 */
	EnhancedNumberValidator(S min, S max, S step):GenericNumberValidator<S>(unrecognizedId, min, max, step){}

	/**
	 * Constructs an EnhancedNumberValidator without explicit minimums or maximums.
	 *
	 * @param step The increments at which the value should be changed. This is mostly used for 
	 * the QSpinBox that is used in the Optika GUI. If you're not using the GUI, you may ignore this parameter.
	 */
	EnhancedNumberValidator(S step):GenericNumberValidator<S>(unrecognizedId, step){}
};

/**
 * A specific validator used to validate entry's of type int.
 */
template <>
class EnhancedNumberValidator<int> : public GenericNumberValidator<int>{
public:
	/**
	 * Construcsts an EnhancedNumberValidator of type int with no
	 * minimum or maximum.
	 */
	EnhancedNumberValidator():GenericNumberValidator<int>(intId, intDefaultStep){}

	/**
	 * Constructs an Enhanced number validator for type int.
	 *
	 * @param min The minimum acceptable value for this validator.
	 * @param max The maximum acceptable value for this validator.
	 * @param step The increments at which the value should be changed. This is mostly used for 
	 * the QSpinBox that is used in the Optika GUI. If you're not using the GUI, you may ignore this parameter.
	 */
	EnhancedNumberValidator(int min, int max, int step=intDefaultStep):GenericNumberValidator<int>(intId, min, max, step){}

	/**
	 * Applies an EnhancedNumberValidator of type int to a QSpinBox
	 *
	 * @param validator The validator to be useed.
	 * @param spinBox The spinBox on which to apply the validator.
	 */
	static void applyToSpinBox(Teuchos::RCP<const EnhancedNumberValidator<int> > validator, QSpinBox *spinBox){
		if(!Teuchos::is_null(validator)){
			spinBox->setMinimum(validator->min());
			spinBox->setMaximum(validator->max());
			spinBox->setSingleStep(validator->getStep());
		}
		else{
			spinBox->setMinimum(INT_MIN);
			spinBox->setMaximum(INT_MAX);
			spinBox->setSingleStep(intDefaultStep);
		}
	}
};

/**
 * A specific validator used to validate values of type short.
 */
template<>
class EnhancedNumberValidator<short> : public GenericNumberValidator<short>{
public:
	/**
	 * Construcsts an EnhancedNumberValidator of type short with no
	 * minimum or maximum.
	 */
	EnhancedNumberValidator():GenericNumberValidator<short>(shortId, shortDefaultStep){}

	/**
	 * Constructs an EnhancedNumberValidator of type short.
	 *
	 * @param min The minimum acceptable value for this validator.
	 * @param max The maximum acceptable value for this validator.
	 * @param step The increments at which the value should be changed. This is mostly used for 
	 * the QSpinBox that is used in the Optika GUI. If you're not using the GUI, you may ignore this parameter.
	 */
	EnhancedNumberValidator(short min, short max, short step=shortDefaultStep):GenericNumberValidator<short>(shortId, min, max, step){}

	/**
	 * Applies an EnhancedNumberValidator of type short to a QSpinBox
	 *
	 * @param validator The validator to be useed.
	 * @param spinBox The spinBox on which to apply the validator.
	 */
	static void applyToSpinBox(Teuchos::RCP<const EnhancedNumberValidator<short> > validator, QSpinBox *spinBox){
		if(!Teuchos::is_null(validator)){
			spinBox->setMinimum(validator->min());
			spinBox->setMaximum(validator->max());
			spinBox->setSingleStep(validator->getStep());
		}
		else{
			spinBox->setMinimum(SHRT_MIN);
			spinBox->setMaximum(SHRT_MAX);
		}
	}
};

/**
 * A specific validator used to validate entry's of type long long int.
 */
/*template<long long int>
class EnhancedNumberValidator<long long int> : public GenericNumberValidator<long long int>{
public:
	**
	 * Construcsts an EnhancedNumberValidator of type long long with no
	 * minimum or maximum.
	 *
	EnhancedNumberValidator():GenericNumberValidator(longlongId, longlongDefaultStep){}
	**
	 * Constructs an Enhanced number validator for type long long int.
	 *
	 * @param min The minimum acceptable value for this validator.
	 * @param max The maximum acceptable value for this validator.
	 * @param step The increments at which the value should be changed. This is mostly used for 
	 * the QSpinBox that is used in the Optika GUI. If you're not using the GUI, you may ignore this parameter.
	 *
	EnhancedNumberValidator(long long int min, long long int max, long long int step=longlongDefaultStep):
	:GenericNumberValidator<long long int>(longlongId, min, max, step){}

	**
	 * Applies an EnhancedNumberValidator of type long long int to a QwwLongSpinBox
	 *
	 * @param validator The validator to be useed.
	 * @param spinBox The spinBox on which to apply the validator.
	 *
	static void applyToSpinBox(Teuchos::RCP<const EnhancedNumberValidator<long long int> > validator, QwwLongSpinBox *spinBox){
		if(!Teuchos::is_null(validator)){
			spinBox->setMinimum(tempMin);
			spinBox->setMaximum(validator->max());
			spinBox->setSingleStep(validator->getStep());
		}
		else{
			spinBox->setMinimum(LLONG_MIN);
			spinBox->setMaximum(LLONG_MAX);
			spinBox->setSingleStep(longlongDefaultStep);
		}
	}
};*/

/**
 * A specific validator used to validate values of type double.
 */
template<>
class EnhancedNumberValidator<double> : public GenericNumberValidator<double>{
public:
	/**
	 * Construcsts an EnhancedNumberValidator of type double with no
	 * minimum or maximum.
	 */
	EnhancedNumberValidator():GenericNumberValidator<double>(doubleId, doubleDefaultStep), precision(doubleDefaultPrecision){}

	/**
	 * Constructs an EnhancedNumberValidator of type double.
	 *
	 * @param min The minimum acceptable value for this validator.
	 * @param max The maximum acceptable value for this validator.
	 * @param step The increments at which the value should be changed. This is mostly used for 
	 * the QSpinBox that is used in the Optika GUI. If you're not using the GUI, you may ignore this parameter.
	 * @param precision This determines the precision at which the number should be displayed in the GUI. 
	 * NOTE: THIS DOES NOT ACTUALLY SPECIFY THE PRECISION USED IN STORING THE VARIABLE. IT IS FOR GUI PURPOSES ONLY!
	 */
	EnhancedNumberValidator(double min, double max, double step=doubleDefaultStep, int precision=doubleDefaultPrecision)
	:GenericNumberValidator<double>(doubleId, min, max, step), precision(precision){}


	/**
	 * Sets the precision specified for the validator. 
	 * NOTE: THIS IS PRECISION USED BY THE GUI. NOT THE ACTUAL PRECISION USED
	 * TO STORE THE VARIABLE.
	 *
	 * @param The precision specific for the validator.
	 */
	void setPrecision(int precision){
		this->precision = precision;
	}

	/**
	 * Gets the precision specified for the validator. 
	 * NOTE: THIS IS PRECISION USED BY THE GUI. NOT THE ACTUAL PRECISION USED
	 * TO STORE THE VARIABLE.
	 *
	 * @return The precision specific for the validator.
	 */
	int getPrecision() const{
		return precision;
	}

	/**
	 * Applies an EnhancedNumberValidator of type double to a QDoubleSpinBox
	 *
	 * @param validator The validator to be useed.
	 * @param spinBox The SpinBox on which to apply the validator.
	 */
	static void applyToSpinBox(Teuchos::RCP<const EnhancedNumberValidator<double> > validator, QDoubleSpinBox *spinBox){
		if(!Teuchos::is_null(validator)){
			spinBox->setMinimum(validator->min());
			spinBox->setMaximum(validator->max());
			spinBox->setSingleStep(validator->getStep());
			spinBox->setDecimals(validator->getPrecision());
		}
		else{
			spinBox->setMinimum(-DBL_MAX);
			spinBox->setMaximum(DBL_MAX);
			spinBox->setSingleStep(doubleDefaultStep);
			spinBox->setDecimals(doubleDefaultPrecision);
		}
	}

private:
	/**
	 * The precision with which the variable should be displayed in a GUI.
	 * NOTE: THIS IS PRECISION USED BY THE GUI. NOT THE ACTUAL PRECISION USED
	 * TO STORE THE VARIABLE.
	 */
	int precision;
};

/**
 * A specific validator used to validat values of type float.
 */
template<>
class EnhancedNumberValidator<float> : public GenericNumberValidator<float>{
public:

	/**
	 * Construcsts an EnhancedNumberValidator of type float with no
	 * minimum or maximum.
	 */
	EnhancedNumberValidator():GenericNumberValidator<float>(floatId, floatDefaultStep), precision(floatDefaultPrecision){}

	/**
	 * Constructs an EnhancedNumberValidator of type float.
	 *
	 * @param min The minimum acceptable value for this validator.
	 * @param max The maximum acceptable value for this validator.
	 * @param step The increments at which the value should be changed. This is mostly used for 
	 * the QSpinBox that is used in the Optika GUI. If you're not using the GUI, you may ignore this parameter.
	 * @param precision This determines the precision at which the number should be displayed in the GUI. 
	 * NOTE: THIS DOES NOT ACTUALLY SPECIFY THE PRECISION USED IN STORING THE VARIABLE. IT IS FOR GUI PURPOSES ONLY!
	 */
	EnhancedNumberValidator(float min, float max, float step=floatDefaultStep, int precision=floatDefaultPrecision)
	:GenericNumberValidator<float>(floatId, min, max, step), precision(precision){}

	/**
	 * Sets the precision specified for the validator. 
	 * NOTE: THIS IS PRECISION USED BY THE GUI. NOT THE ACTUAL PRECISION USED
	 * TO STORE THE VARIABLE.
	 *
	 * @param The precision specific for the validator.
	 */
	void setPrecision(int precision){
		this->precision = precision;
	}

	/**
	 * Gets the precision specified for the validator. 
	 * NOTE: THIS IS PRECISION USED BY THE GUI. NOT THE ACTUAL PRECISION USED
	 * TO STORE THE VARIABLE.
	 *
	 * @return The precision specific for the validator.
	 */
	int getPrecision() const{
		return precision;
	}

	/**
	 * Applies an EnhancedNumberValidator of type float to a QDoubleSpinBox.
	 *
	 * @param validator The validator to be useed.
	 * @param spinBox The SpinBox on which to apply the validator.
	 */
	static void applyToSpinBox(Teuchos::RCP<const EnhancedNumberValidator<float> > validator, QDoubleSpinBox *spinBox){
		if(!Teuchos::is_null(validator)){
			spinBox->setMinimum(validator->min());
			spinBox->setMaximum(validator->max());
			spinBox->setSingleStep(validator->getStep());
			spinBox->setDecimals(validator->getPrecision());
		}
		else{
			spinBox->setMinimum(-FLT_MAX);
			spinBox->setMaximum(FLT_MAX);
			spinBox->setSingleStep(floatDefaultStep);
			spinBox->setDecimals(floatDefaultPrecision);
		}
	}

private:
	/**
	 * The precision with which the variable should be displayed in a GUI.
	 * NOTE: THIS IS PRECISION USED BY THE GUI. NOT THE ACTUAL PRECISION USED
	 * TO STORE THE VARIABLE.
	 */
	int precision;
}; 

/**
 * Simply indicates that the parameter entry with this validator should
 * contain a filename.
 */
class FileNameValidator : public Teuchos::ParameterEntryValidator{
public:
	/**
	 * Constructs a FileNameValidator.
	 *
	 * @param mustAlreadyExist True if the file the user specifies should already exists, false otherwise.
	 */
	FileNameValidator(bool mustAlreadyExist=false);

	/**
	 * Gets the variable describing whether or not this validator wants the file that is specified to
	 * already exist.
	 *
	 * @return Whether or not the validator requires the file to already exist
	 */
	bool fileMustExist() const;

	/**
	 * Sets whether or not the validator requires the file to already exist.
	 *
	 * @param shouldFileExist True if the file should already exist, false otherwise.
	 * @return The new value of the shouldFileExist variable.
	 */
	bool setFileMustExist(bool shouldFileExist);

	Teuchos::RCP<const Teuchos::Array<std::string> > validStringValues() const;

	void validate(Teuchos::ParameterEntry const &entry, std::string const &paramName, std::string const &sublistName) const;

	void printDoc(std::string const &docString, std::ostream &out) const;
private:
	/**
	 * Whether or not the file specified in the parameter should already exist.
	 */
	bool mustAlreadyExist;
};

/**
 * A simple validator that only allows certain string values to be choosen.
 */
class StringValidator : public Teuchos::ParameterEntryValidator{
public:
	typedef Teuchos::Array<std::string> ValueList;
	/**
	 * Constructs a StringValidator.
	 */
	StringValidator(ValueList validStrings);

	/**
	 * Sets the Array of valid strings and returns what the current array of valid
	 * string now is.
	 *
	 * @param validStrings What the array for the valid strings should contain.
	 * @return What the arry for the valid strings now conatians.
	 */
	const ValueList setValidStrings(ValueList validStrings);

	Teuchos::RCP<const Teuchos::Array<std::string> > validStringValues() const;

	void validate(Teuchos::ParameterEntry const &entry, std::string const &paramName, std::string const &sublistName) const;

	void printDoc(std::string const &docString, std::ostream &out) const;
private:
	/**
	 * An array containing a list of all the valid string values.
	 */
	ValueList validStrings;
};

/**
 * An Abstract base class for all ArrayValidators
 */
class ArrayValidator : public Teuchos::ParameterEntryValidator{
public:
	/**
	 * Constructs a ArrayValidator.
	 *
	 * @param prototypeValidator The validator to be used on each
	 * entry in the array.
	 */
	ArrayValidator(Teuchos::RCP<Teuchos::ParameterEntryValidator> prototypeValidator):ParameterEntryValidator(),
		       prototypeValidator(prototypeValidator){}

	virtual Teuchos::RCP<const Teuchos::Array<std::string> > validStringValues() const{
		return prototypeValidator->validStringValues();
	}
	virtual void validate(Teuchos::ParameterEntry const &entry, std::string const &paramName, std::string const &sublistName) const =0;
	virtual void printDoc(std::string const &docString, std::ostream &out) const =0;
protected:
	/**
	 * The prototype validator to be applied to each entry in the Array.
	 */
	Teuchos::RCP<Teuchos::ParameterEntryValidator> prototypeValidator;
};

/**
 * A wrapper class for Arrays using a StringToIntegralParameterEntryValidator.
 */
class ArrayStringValidator: public ArrayValidator{
public:
	/**
	 * Constructs an ArrayStringValidator
	 *
	 * @param prototypeValidator The ParameterEntry validator containing a list of valid string values to be used 
	 * on each entry in the Array.
	 */
	ArrayStringValidator(Teuchos::RCP<Teuchos::ParameterEntryValidator> prototypeValidator)
	:ArrayValidator(prototypeValidator){}
	
	/**
	 * Retruns the prototype validator being used by the ArrayStringValidator.
	 *
	 * @return The prototype validator being used by the ArrayStringValidator.
	 */
	Teuchos::RCP<const Teuchos::ParameterEntryValidator> getPrototype() const{
		return prototypeValidator;
	}

	void validate(Teuchos::ParameterEntry const &entry, std::string const &paramName, std::string const &sublistName) const{
		Teuchos::any anyValue = entry.getAny(true);
		if(anyValue.type() == typeid(Teuchos::Array<std::string>)){
			Teuchos::Array<std::string> extracted = Teuchos::getValue<Teuchos::Array<std::string> >(entry);
			
			std::string currentString;
			Teuchos::RCP< const Teuchos::Array<std::string> > validStrings = validStringValues();
			for(int i = 0; i<extracted.size(); ++i){
				currentString = extracted[i];
				Teuchos::Array<std::string>::const_iterator it = std::find(validStrings->begin(), validStrings->end(), currentString);
				if(it == validStrings->end()){
					std::stringstream oss;
					std::string msg;
					oss << "Aww shoot! Sorry bud, but it looks like the \"" << paramName << "\"" <<
					" parameter in the \"" << sublistName << "\" sublist didn't quite work out.\n" <<
					"No need to fret though. I'm sure it's just a small mistake. Maybe the information below "<<
					"can help you figure out what went wrong.\n\n"
					"Error: The value that was entered at " << i << " in the array does't fall within " <<
					"the rang set by the validtor.\n" <<
					"Parameter: " << paramName << "\n" << 
					"Value entered at " << i << ": "<<
					extracted[i] << "\n" <<
					"Exceptable Values:\n";
					for(int j=0; j<validStrings->size(); ++j){
						oss << "	" << (*validStrings)[j] << "\n";
					}
					msg = oss.str();
					throw Teuchos::Exceptions::InvalidParameterValue(msg);
				}	
			}
		}
		else{
			const std::string &entryName = entry.getAny(false).typeName();
			std::stringstream oss;
			std::string msg;
			oss << "Aww shoot! Sorry bud, but it looks like the \"" << paramName << "\"" <<
			" parameter in the \"" << sublistName << "\" sublist didn't quite work out.\n" <<
			"No need to fret though. I'm sure it's just a small mistake. Maybe the information below "<<
			"can help you figure out what went wrong.\n\n"
			"Error: The value you entered was the wrong type.\n" <<
			"Parameter: " << paramName << "\n" <<
			"Type specified: " << entryName << "\n" <<
			"Type accepted: " << typeid(Teuchos::Array<std::string>).name() << "\n";
			msg = oss.str();
			throw Teuchos::Exceptions::InvalidParameterType(msg);
		}

	}

	void printDoc(std::string const &docString, std::ostream &out) const{
		Teuchos::StrUtils::printLines(out,"# ",docString);
		std::string toPrint;
		toPrint += "Prototype Validator:\n";
		prototypeValidator->printDoc(toPrint, out);
	}
};


/**
 * A wrapper class allowing EnhancedNumberValidators to be applied to arrays.
 */
template <class S>
class ArrayNumberValidator: public ArrayValidator{
public:
	/**
	 * Constructs and ArrayNumberValidator.
	 *
	 * @param prototypeValidator The validator to be applied to each entry
	 * in the array.
	 */
	ArrayNumberValidator(Teuchos::RCP<EnhancedNumberValidator<S> > prototypeValidator)
	:ArrayValidator(prototypeValidator){}

	/**
	 * Retruns the prototype validator being used by the ArrayNumberValidator.
	 *
	 * @return The prototype validator being used by the ArrayNumberValidator.
	 */
	Teuchos::RCP<const EnhancedNumberValidator<S> > getPrototype() const{
		return Teuchos::rcp_static_cast<EnhancedNumberValidator<S> > (prototypeValidator);
	}

	void validate(Teuchos::ParameterEntry const &entry, std::string const &paramName, std::string const &sublistName) const{
		Teuchos::any anyValue = entry.getAny(true);
		if(anyValue.type() == typeid(Teuchos::Array<S>)){
			Teuchos::Array<S> extracted = Teuchos::any_cast<Teuchos::Array<S> >(anyValue);
			for(int i = 0; i<extracted.size(); ++i){
				if(!( extracted[i] >= getPrototype()->min() &&  extracted[i] <= getPrototype()->max())){
					std::stringstream oss;
					std::string msg;
					oss << "Aww shoot! Sorry bud, but it looks like the \"" << paramName << "\"" <<
					" parameter in the \"" << sublistName << "\" sublist didn't quite work out.\n" <<
					"No need to fret though. I'm sure it's just a small mistake. Maybe the information below "<<
					"can help you figure out what went wrong.\n\n"
					"Error: The value that was entered at \"" << i << "\" in the array does't fall within " <<
					"the rang set by the validtor.\n" <<
					"Parameter: " << paramName << "\n" <<
					"Min: " << getPrototype()->min() << "\n" <<
					"Max: " << getPrototype()->max() << "\n" <<
					"Value entered at " << i << ": "<<
					extracted[i] << "\n";
					msg = oss.str();
					throw Teuchos::Exceptions::InvalidParameterValue(msg);
				}	
			}
		}
		else{
			const std::string &entryName = entry.getAny(false).typeName();
			std::stringstream oss;
			std::string msg;
			oss << "Aww shoot! Sorry bud, but it looks like the \"" << paramName << "\"" <<
			" parameter in the \"" << sublistName << "\" sublist didn't quite work out.\n" <<
			"No need to fret though. I'm sure it's just a small mistake. Maybe the information below "<<
			"can help you figure out what went wrong.\n\n"
			"Error: The value you entered was the wrong type.\n" <<
			"Parameter: " << paramName << "\n" <<
			"Type specified: " << entryName << "\n" <<
			"Type accepted: " << typeid(Teuchos::Array<S>).name() << "\n";
			msg = oss.str();
			throw Teuchos::Exceptions::InvalidParameterType(msg);
		}

	}

	void printDoc(std::string const &docString, std::ostream &out) const{
		Teuchos::StrUtils::printLines(out,"# ",docString);
		std::string toPrint;
		toPrint += "Prototype Validator:\n";
		prototypeValidator->printDoc(toPrint, out);
	}
};


/**
 * A wrapper class allowing a FileNameValidator to be applied to an array.
 */
class ArrayFileNameValidator : public ArrayValidator{
public:
	/**
	 * Constructs a FileNameArrayValidator.
	 *
	 * @param prototypeValidator The validator to be applied to each entry in the array.
	 */
	ArrayFileNameValidator(Teuchos::RCP<FileNameValidator> prototypeValidator)
	:ArrayValidator(prototypeValidator){}

	/**
	 * Retruns the prototype validator being used by the ArrayFileNameValidator.
	 *
	 * @return The prototype validator being used by the ArrayFileNameValidator.
	 */
	Teuchos::RCP<const FileNameValidator> getPrototype() const{
		return Teuchos::rcp_static_cast<FileNameValidator>(prototypeValidator);
	}

	void validate(Teuchos::ParameterEntry const &entry, std::string const &paramName, std::string const &sublistName) const{
		Teuchos::any anyValue = entry.getAny(true);
		if(!(anyValue.type() == typeid(Teuchos::Array<std::string>))){
			const std::string &entryName = entry.getAny(false).typeName();
			std::stringstream oss;
			std::string msg;
			oss << "Aww shoot! Sorry bud, but it looks like the \"" << paramName << "\"" <<
			" parameter in the \"" << sublistName << "\" sublist didn't quite work out.\n" <<
			"No need to fret though. I'm sure it's just a small mistake. Maybe the information below "<<
			"can help you figure out what went wrong.\n\n"
			"Error: The value you entered was the wrong type.\n" <<
			"Parameter: " << paramName << "\n" <<
			"Type specified: " << entryName << "\n" <<
			"Type accepted: " << anyValue.typeName() << "\n";
			msg = oss.str();
			throw Teuchos::Exceptions::InvalidParameterType(msg);
		}
		else if(getPrototype()->fileMustExist()){
			Teuchos::Array<std::string> extracted = Teuchos::any_cast<Teuchos::Array<std::string> >(anyValue);
			for(int i = 0; i<extracted.size(); ++i){
				std::string fileName = extracted[i];
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
					"File name specified index " << i <<": " << fileName << "\n";
					msg = oss.str();
					throw Teuchos::Exceptions::InvalidParameterValue(msg);
				}
			}
		}

	}

	void printDoc(std::string const &docString, std::ostream &out) const{
		Teuchos::StrUtils::printLines(out,"# ",docString);
		std::string toPrint;
		toPrint += "Prototype Validator:\n";
		prototypeValidator->printDoc(toPrint, out);
	}
};





}
#endif // OPTIKA_SPECIFICPARAMETERENTRYVALIDATORS_HPP_
