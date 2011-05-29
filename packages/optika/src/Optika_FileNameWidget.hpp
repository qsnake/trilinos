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
#ifndef OPTIKA_FILENAMEWIDGET_HPP_
#define OPTIKA_FILENAMEWIDGET_HPP_
#include <QWidget>
class QLabel;

namespace Optika{

/**
 * A small widget consisting of a button and label that allows the user
 *  to select a file through a QFileDialog. The label displays
 *  the currently selected file.
 */
class FileNameWidget : public QWidget{
	Q_OBJECT
public:
	/**
	 * Constructs a FileNameWidget
	 *
	 * @param currentFileName The Filename with which the widget should be 
	 * initialized.
	 * @param parent The parent widget.
	 */
	FileNameWidget(QString currentFileName=QString(), bool mustAlreadyExist=false, QWidget* parent=0);

	/**
	 * Gets the current filename in the widget.
	 *
	 * @return The current filename in the widget.
	 */
	QString getCurrentFileName();

	/**
	 * Sets the current filename in the widget.
	 *
	 * @param The name to which the widget should be set.
	 */
	void setCurrentFileName(QString newName);

public slots:
	/**
	 * Opens a QFileDialog allowing the user to select a new filename.
	 */
	void getNewFileName();

private:
	/**
	 * The current file name stored in the list.
	 */
	QString currentFileName;

	/**
	 * The label describing the file path.
	 */
	QLabel *pathLabel;
	
	/**
	 * Whether or not the file name specified must already exist.
	 */
	bool mustAlreadyExist;
};

}

#endif //OPTIKA_FILENAMEWIDGET_HPP_
