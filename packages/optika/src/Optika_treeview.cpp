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
#include "Optika_treeview.hpp"
#include <QMessageBox>
namespace Optika{


TreeView::TreeView(TreeModel *treeModel, Delegate *delegate):QTreeView(){
	setModel(treeModel);
	setItemDelegateForColumn(1, delegate);
	if(treeModel->hasDependencies()){
		connect(treeModel, SIGNAL(hideData(int, const QModelIndex&)), this, SLOT(hideRow(int, const QModelIndex&)));
		connect(treeModel, SIGNAL(showData(int, const QModelIndex&)), this, SLOT(showRow(int, const QModelIndex&)));
		connect(treeModel, SIGNAL(badValue(QModelIndex, QString)), this, SLOT(handleBadValue(QModelIndex, QString)));
		connect(delegate, SIGNAL(closeEditor(QWidget*, QAbstractItemDelegate::EndEditHint)), this, SLOT(checkForOtherBadValues()));
		treeModel->issueInitilizationSignals();
	}
	setAnimated(true);
	setAlternatingRowColors(true);
}

void TreeView::showRow(int row, const QModelIndex& parent){
	if(isRowHidden(row, parent)){
		setRowHidden(row, parent, false);
	}
}

void TreeView::hideRow(int row, const QModelIndex& parent){
	if(!isRowHidden(row, parent)){
		setRowHidden(row, parent, true);
	}
}

void TreeView::handleBadValue(QModelIndex badValueIndex, QString message){
	if(state() != EditingState && !isRowHidden(badValueIndex.row(), badValueIndex.parent())){
		QMessageBox::warning(this, "Bad parameter value", message);
		setCurrentIndex(badValueIndex);
		edit(badValueIndex);
	}
	else if(!isRowHidden(badValueIndex.row(), badValueIndex.parent())){
		invalidInicies.enqueue(invalidIndex(badValueIndex, message));
	}
}

void TreeView::checkForOtherBadValues(){
	if(invalidInicies.size() != 0){
		invalidIndex needsToBeEdited = invalidInicies.dequeue();
		handleBadValue(needsToBeEdited.first, needsToBeEdited.second);
	}
}

}

