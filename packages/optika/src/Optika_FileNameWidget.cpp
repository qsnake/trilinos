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
#include "Optika_FileNameWidget.hpp"
#include <QPushButton>
#include <QLabel>
#include <QVBoxLayout>
#include <QFileDialog>


namespace Optika{


FileNameWidget::FileNameWidget(QString currentFileName, bool mustAlreadyExist, QWidget *parent)
	:QWidget(parent),
	currentFileName(currentFileName),
	mustAlreadyExist(mustAlreadyExist)
{
	QPushButton *changeButton = new QPushButton("Change Path",this);
	connect(changeButton, SIGNAL(clicked(bool)), this, SLOT(getNewFileName()));
	pathLabel = new QLabel(currentFileName,this);
	QVBoxLayout *layout = new QVBoxLayout(this);
	layout->addWidget(changeButton);
	layout->addWidget(pathLabel);
	setLayout(layout);
}

QString FileNameWidget::getCurrentFileName(){
	return currentFileName;
}

void FileNameWidget::setCurrentFileName(QString newName){
	currentFileName = newName;
	pathLabel->setText(newName);	
}

void FileNameWidget::getNewFileName(){
	QString defaultPath;
	if(currentFileName == ""){
		defaultPath = QDir::homePath();
	}
	else{
		defaultPath = currentFileName;
	}
	if(mustAlreadyExist){
		setCurrentFileName(QFileDialog::getOpenFileName(this, tr("File"), defaultPath));
	}
	else{
		setCurrentFileName(QFileDialog::getSaveFileName(this, tr("File"), defaultPath));
	}
}


}

