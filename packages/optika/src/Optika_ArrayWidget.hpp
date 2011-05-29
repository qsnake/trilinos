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
#ifndef OPTIKA_ARRAYWIDGET_HPP_
#define OPTIKA_ARRAYWIDGET_HPP_

#include <QDialog>
#include <QModelIndex>
#include <QPushButton>
#include <QGridLayout>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QLineEdit>
#include <QGridLayout>
#include <QScrollArea>
#include <QLabel>
#include <vector>
#include "Optika_treemodel.hpp"
#include "Optika_FileNameWidget.hpp"

namespace Optika {


/**
 * A templated abstract base class for all other array editing widgets. Note the absence of the Q_OBJECT
 * macro. This is becuase classes using the Q_OBJECT macro can't be templated (bummer). The macro is therfore
 * present in the subclasses.
 */ 
template <class S>
class GenericArrayWidget : public QDialog{
public:
	/**
	 * Constructs a GenericArrayWidget.
	 *
	 * @param index The index of the array that is being edited.
	 * @param type The type of the array.
	 * @param parent The parent widget.
	 */
	GenericArrayWidget(const QModelIndex index, QString type, QWidget *parent=0):
		QDialog(parent),
		model((TreeModel*)index.model()),
		index(index), 
		entryValidator(model->getValidator(index)),
		baseArray(model->getArray<S>(index)),
		type(type)
	{
		setModal(true);
		setSizeGripEnabled(true);
		arrayContainer = new QWidget(this);

		QScrollArea *scrollArea = new QScrollArea(this);
		scrollArea->setWidget(arrayContainer);
		scrollArea->setWidgetResizable(true);

		QPushButton *doneButton = new QPushButton(tr("Done"));
		QPushButton *cancelButton = new QPushButton(tr("Cancel"));
		connect(doneButton, SIGNAL(clicked(bool)), this, SLOT(accept()));
		connect(cancelButton, SIGNAL(clicked(bool)), this, SLOT(reject()));
		QGridLayout *layout = new QGridLayout(this);
		layout->addWidget(scrollArea,0,0,1,3);
		layout->addWidget(doneButton,1,2);
		layout->addWidget(cancelButton,1,1);

		this->setLayout(layout);

		setWindowTitle(model->data(index.sibling(index.row(),0),Qt::DisplayRole).toString());
	}
	
	/**
	 * Gets the type of array being edited.
	 *
	 * @return The type of array being edited.
	 */
	QString getType(){
		return type;
	}

	/**
	 * Gets a string representing what should be saved back to the model.
	 * When reimplemented in a subclass, it should be a slot.
	 *
	 * @return A string representing what should be saved back to the model.
	 */
	virtual std::string saveData() = 0;

	/**
	 * Sets all of the values in the array widget to what they initially should be.
	 * When reimplemented in a subclass, it should be a slot.
	 *
	 * @param values The values to which the array should be set.
	 */
	virtual void initializeValues(QString values) = 0;

	/**
	 * Called when the user has entered in their desired values and is done editing
	 * the array. When reimplemented in a subclass, it should be a slot.
	 */
	virtual void accept() =0;

protected:
	/**
	 * Convienece typedef. Represents an array of QWidgets.
	 */
	typedef std::vector<QWidget*> WVector;

	/**
	 * Conatins the editing widgets (e.g. QLineEdits and QSpinBoxes) comprising the array editor.
	 */
	WVector widgetVector;

	/**
	 * A pointer to the TreeModel the ArrayWidget is editing.
	 */
	TreeModel *model;

	/**
	 * The index in the TreeModel of the parameter the ArrayWidget is
	 * editing.
	 */
	QModelIndex index;

	/**
	 * The validator being used on the array.
	 */
	Teuchos::RCP<const Teuchos::ParameterEntryValidator> entryValidator;	

	/**
	 * Sets up the layout for the arrayContainer, including adding what ever editing
	 * widget should be used for the particual type of array.
	 */
	void setupArrayLayout(){
		QGridLayout *widgetLayout = new QGridLayout(arrayContainer);
		for(int i=0; i<baseArray.size(); ++i){
			widgetLayout->addWidget(new QLabel("Item: " +QString::number(i)),0,i,Qt::AlignLeft);
			QWidget* editorWidget = getEditorWidget();
			widgetLayout->addWidget(editorWidget,1,i,Qt::AlignLeft);
			widgetVector.push_back(editorWidget);
		}
		arrayContainer->setLayout(widgetLayout);
	}

private:
	/**
	 * The widget containing all of the editing widgets (e.g.
	 * QLineEdits, and QSpinBoxes) that comprise the array editor.
	 */
	QWidget *arrayContainer;
	
	/**
	 * The array to be edited.
	 */
	Teuchos::Array<S> baseArray;

	/**
	 * The type of array.
	 */
	QString type;
	
	/**
	 * Gets the widget to be used as an editor for each entry in the array.
	 */
	virtual QWidget* getEditorWidget() = 0;
};

/**
 * A widget for editing Arrays of type int.
 */
class IntArrayWidget: public GenericArrayWidget<int>{
	Q_OBJECT
public:
	/**
	 * Constructs an IntArrayWidget.
	 *
	 * @param index The index of the array that is being edited.
	 * being edited.
	 * @param type The type of the array.
	 * @param parent The parent widget.
	 */
	IntArrayWidget(const QModelIndex index, QString type, QWidget *parent=0):GenericArrayWidget<int>(index, type, parent){
		setupArrayLayout();
		initializeValues(index.model()->data(index).toString()); 
	}

public slots:
	void accept(){
		model->setData(index, QString::fromStdString(saveData()), Qt::EditRole);
		done(QDialog::Accepted);
	}

	std::string saveData(){
		Teuchos::Array<int> toReturn;
		for(WVector::iterator it = widgetVector.begin(); it != widgetVector.end(); ++it){
			toReturn.push_back(((QSpinBox*)(*it))->value());
		}
		return toReturn.toString();
	}

	void initializeValues(QString values){
		QStringList valueList = getValues(values); 
		int i =0;
		for(WVector::iterator it = widgetVector.begin(); it != widgetVector.end(); ++it, ++i){
			static_cast<QSpinBox*>(*it)->setValue(valueList.at(i).toInt());
		}

	}

private:
	QWidget* getEditorWidget(){
		QSpinBox *newSpin = new QSpinBox(this);
		Teuchos::RCP<const EnhancedNumberValidator<int> > validator = Teuchos::null;
		if(!Teuchos::is_null(entryValidator)){
			validator = Teuchos::rcp_dynamic_cast<const ArrayNumberValidator<int> >(entryValidator,true)->getPrototype();
		}
		EnhancedNumberValidator<int>::applyToSpinBox(validator, newSpin);
		return newSpin;
	}
};

/**
 * A widget for editing Arrays of type short.
 */
class ShortArrayWidget: public GenericArrayWidget<short>
{
	Q_OBJECT
public:
	/**
	 * Constructs a ShortArrayWidget.
	 *
	 * @param index The index of the array that is being edited.
	 * being edited.
	 * @param type The type of the array.
	 * @param parent The parent widget.
	 */
	ShortArrayWidget(const QModelIndex index, QString type, QWidget *parent=0):GenericArrayWidget<short>(index, type, parent){
		setupArrayLayout();
		initializeValues(index.model()->data(index).toString()); 
	}

public slots:
	void accept(){
		model->setData(index, QString::fromStdString(saveData()), Qt::EditRole);
		done(QDialog::Accepted);
	}

	std::string saveData(){
		Teuchos::Array<short> toReturn;
		for(WVector::iterator it = widgetVector.begin(); it != widgetVector.end(); ++it){
			toReturn.push_back(((QSpinBox*)(*it))->value());
		}
		return toReturn.toString();
	}

	void initializeValues(QString values){
		QStringList valueList = getValues(values); 
		int i =0;
		for(WVector::iterator it = widgetVector.begin(); it != widgetVector.end(); ++it, ++i){
			static_cast<QSpinBox*>(*it)->setValue((short)valueList.at(i).toInt());
		}

	}

private:
	QWidget* getEditorWidget(){
		QSpinBox *newSpin = new QSpinBox(this);
		Teuchos::RCP<const EnhancedNumberValidator<short> > validator = Teuchos::null;
		if(!Teuchos::is_null(entryValidator)){
			validator = Teuchos::rcp_dynamic_cast<const ArrayNumberValidator<short> >(entryValidator,true)->getPrototype();
		}
		EnhancedNumberValidator<short>::applyToSpinBox(validator, newSpin);
		return newSpin;
	}
};

/**
 * A widget for editing Arrays of type long long int.
 */
/*
class LongLongArrayWidget: public GenericArrayWidget<long long int>{
	Q_OBJECT
public:
	**
	 * Constructs an LongLongArrayWidget.
	 *
	 * @param index The index of the array that is being edited.
	 * @param type The type of the array.
	 * @param parent The parent widget.
	 *
	LongLongArrayWidget(const QModelIndex index, QString type, QWidget *parent=0):GenericArrayWidget<long long int>(index, type, parent){
		setupArrayLayout();
		initializeValues(index.model()->data(index).toString()); 
	}

public slots:
	void accept(){
		model->setData(index, QString::fromStdString(saveData()), Qt::EditRole);
		done(QDialog::Accepted);
	}

	std::string saveData(){
		Teuchos::Array<long long int> toReturn;
		for(WVector::iterator it = widgetVector.begin(); it != widgetVector.end(); ++it){
			toReturn.push_back(((QwwLongSpinBox*)(*it))->value());
		}
		return toReturn.toString();
	}

	void initializeValues(QString values){
		QStringList valueList = getValues(values); 
		int i =0;
		for(WVector::iterator it = widgetVector.begin(); it != widgetVector.end(); ++it, ++i){
			static_cast<QwwLongSpinBox*>(*it)->setValue(valueList.at(i).toLongLong());
		}

	}

private:
	QWidget* getEditorWidget(){
		QSpinBox *newSpin = new QwwLongSpinBox(this);
		Teuchos::RCP<const EnhancedNumberValidator<long long int> > validator = null;
		if(!Teuchos::is_null(entryValidator)){
			validator = Teuchos::rcp_dynamic_cast<const EnhancedNumberValidator<long long int> >(entryValidator,true);
		}
		EnhancedNumberValidator<long long int>::applyToSpinBox(validator, newSpin);
		return newSpin;
	}
};*/

/**
 * A widget for editing Arrays of type double.
 */
class DoubleArrayWidget: public GenericArrayWidget<double>
{
	Q_OBJECT
public:
	/**
	 * Constructs a DoubleArrayWidget.
	 *
	 * @param index The index of the array that is being edited.
	 * @param type The type of the array.
	 * @param parent The parent widget.
	 */
	DoubleArrayWidget(const QModelIndex index, QString type, QWidget *parent=0):GenericArrayWidget<double>(index, type, parent){
		setupArrayLayout();
		initializeValues(index.model()->data(index).toString()); 
	}

public slots:
	void accept(){
		model->setData(index, QString::fromStdString(saveData()), Qt::EditRole);
		done(QDialog::Accepted);
	}

	std::string saveData(){
		Teuchos::Array<double> toReturn;
		for(WVector::iterator it = widgetVector.begin(); it != widgetVector.end(); ++it){
			toReturn.push_back(((QDoubleSpinBox*)(*it))->value());
		}
		return toReturn.toString();
	}

	void initializeValues(QString values){
		QStringList valueList = getValues(values); 
		int i =0;
		for(WVector::iterator it = widgetVector.begin(); it != widgetVector.end(); ++it, ++i){
			static_cast<QDoubleSpinBox*>(*it)->setValue(valueList.at(i).toDouble());
		}
	}

private:
	QWidget* getEditorWidget(){
		QDoubleSpinBox *newSpin = new QDoubleSpinBox(this);
		Teuchos::RCP<const EnhancedNumberValidator<double> > validator = Teuchos::null;
		if(!Teuchos::is_null(entryValidator)){
			validator = Teuchos::rcp_dynamic_cast<const ArrayNumberValidator<double> >(entryValidator,true)->getPrototype();
		}
		EnhancedNumberValidator<double>::applyToSpinBox(validator, newSpin);
		return newSpin;
	}
};

/**
 * A widget for editing Arrays of type short.
 */
class FloatArrayWidget: public GenericArrayWidget<float>
{
	Q_OBJECT
public:
	/**
	 * Constructs a FloatArrayWidget.
	 *
	 * @param index The index of the array that is being edited.
	 * @param type The type of the array.
	 * @param parent The parent widget.
	 */
	FloatArrayWidget(const QModelIndex index, QString type, QWidget *parent=0):GenericArrayWidget<float>(index, type, parent){
		setupArrayLayout();
		initializeValues(index.model()->data(index).toString()); 
	}

public slots:
	void accept(){
		model->setData(index, QString::fromStdString(saveData()), Qt::EditRole);
		done(QDialog::Accepted);
	}

	std::string saveData(){
		Teuchos::Array<float> toReturn;
		for(WVector::iterator it = widgetVector.begin(); it != widgetVector.end(); ++it){
			toReturn.push_back(((QDoubleSpinBox*)(*it))->value());
		}
		return toReturn.toString();
	}	

	void initializeValues(QString values){
		QStringList valueList = getValues(values); 
		int i =0;
		for(WVector::iterator it = widgetVector.begin(); it != widgetVector.end(); ++it, ++i){
			static_cast<QDoubleSpinBox*>(*it)->setValue(valueList.at(i).toDouble());
		}

	}

private:
	QWidget* getEditorWidget(){
		QDoubleSpinBox *newSpin = new QDoubleSpinBox(this);
		Teuchos::RCP<const EnhancedNumberValidator<float> > validator = Teuchos::null;
		if(!Teuchos::is_null(entryValidator)){
			validator = Teuchos::rcp_dynamic_cast<const ArrayNumberValidator<float> >(entryValidator,true)->getPrototype();
		}
		EnhancedNumberValidator<float>::applyToSpinBox(validator, newSpin);
		return newSpin;
	}
};

/**
 * A widget for editing an array of strings
 */
class StringArrayWidget : public GenericArrayWidget<std::string>{
	Q_OBJECT
public:
	/**
	 * Constructs a StringArrayWidget.
	 *
	 * @param index The index of the array that is being edited.
	 * @param type The type of the array.
	 * @param parent The parent widget.
	 */
	StringArrayWidget(const QModelIndex index, QString type, QWidget *parent=0):GenericArrayWidget<std::string>(index, type, parent){
		setupArrayLayout();
		initializeValues(index.model()->data(index).toString()); 
	}

	void accept(){
		model->setData(index, QString::fromStdString(saveData()), Qt::EditRole);
		done(QDialog::Accepted);
	}

	std::string saveData(){
		Teuchos::Array<std::string> toReturn;
		for(WVector::iterator it = widgetVector.begin(); it != widgetVector.end(); ++it){
			if(Teuchos::is_null(entryValidator)){
				toReturn.push_back(((QLineEdit*)(*it))->text().toStdString());
			}
			else if(!Teuchos::is_null(Teuchos::rcp_dynamic_cast<const ArrayFileNameValidator>(entryValidator))){
				toReturn.push_back(((FileNameWidget*)(*it))->getCurrentFileName().toStdString());
			}
			else if(entryValidator->validStringValues()->size() !=0){
				toReturn.push_back(((QComboBox*)(*it))->currentText().toStdString());
			}
			else{
				toReturn.push_back(((QLineEdit*)(*it))->text().toStdString());
			}
		}
		return toReturn.toString();
	}	

	void initializeValues(QString values){
		QStringList valueList = getValues(values); 
		int i =0;
		for(WVector::iterator it = widgetVector.begin(); it != widgetVector.end(); ++it, ++i){
			if(Teuchos::is_null(entryValidator)){
				static_cast<QLineEdit*>(*it)->setText(valueList.at(i));
			}
			else if(!Teuchos::is_null(Teuchos::rcp_dynamic_cast<const ArrayFileNameValidator>(entryValidator))){
				static_cast<FileNameWidget*>(*it)->setCurrentFileName(valueList.at(i));
			}
			else if(entryValidator->validStringValues()->size() !=0){
				int currentIndex = static_cast<QComboBox*>(*it)->findText(valueList.at(i));
				if(currentIndex >= 0){
					static_cast<QComboBox*>(*it)->setCurrentIndex(static_cast<QComboBox*>(*it)->findText(valueList.at(i)));
				}
			}
			else{
				static_cast<QLineEdit*>(*it)->setText(valueList.at(i));
			}
		}

	}

private:
	QWidget* getEditorWidget(){
		if(Teuchos::is_null(entryValidator)){
			return new QLineEdit(this);
		}
		else if(!Teuchos::is_null(Teuchos::rcp_dynamic_cast<const ArrayFileNameValidator>(entryValidator))){
			return new FileNameWidget("", Teuchos::rcp_dynamic_cast<const ArrayFileNameValidator>(entryValidator)->getPrototype()->fileMustExist(), this);
		}
		else if(entryValidator->validStringValues()->size() != 0){
			Teuchos::RCP<const Teuchos::Array<std::string> > options = entryValidator->validStringValues();
			QComboBox *newCombo = new QComboBox(this);
			for(Teuchos::Array<std::string>::const_iterator itr = options->begin(); itr != options->end(); ++itr){
				newCombo->addItem(QString::fromStdString(*itr));
			}
			return newCombo;
		}
		else{
			return new QLineEdit(this);
		}
	}
};


}

#endif //OPTIKA_ARRAYWIDGET_HPP_
