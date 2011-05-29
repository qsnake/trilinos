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
#ifndef OPTIKA_TREEVIEW_HPP_
#define OPTIKA_TREEVIEW_HPP_
#include <QTreeView>
#include <QQueue>
#include "Optika_delegate.hpp"
namespace Optika{

class Delegate;
class TreeModel;

/**
 * Class used to view TreeModels
 */
class TreeView : public QTreeView{
	Q_OBJECT
public:
	/**
	 * A pair representing an invalidIndex and why it's invalid
	 */
	typedef std::pair<QModelIndex, QString> invalidIndex;

	/**
	 * Constructs a TreeView.
	 * 
	 * @param treeModel The Tree Model being used with the TreeView.
	 * @param delegate The delegate to be used with the TreeView.
	 */
	TreeView(TreeModel *treeModel, Delegate *delegate);

public slots:
	/**
	 * Used to change the visiblity of a row from hidden to shown.
	 *
	 * @param row The row to be shwon.
	 * @param parent The parent of the item to be shown.
	 */
	void showRow(int row, const QModelIndex& parent);

	/**
	 * Used to change the visiblity of a row from shown to hidden.
	 *
	 * @param row The row to be shwon.
	 * @param parent The parent of the item to be hidden.
	 */
	void hideRow(int row, const QModelIndex& parent);

	/**
	 * Handles any badValue signals that might be emitted by the
	 * TreeModel.
	 *
	 * @param badValueIndex The index of the item with the bad value.
	 * @param A brief message explaining what happened to cause the
	 * treeitem to have an invalid value.
	 */
	void handleBadValue(QModelIndex badValueIndex, QString message);

	/**
	 * Checks to see if there are any other invalid indicies.
	 * If there are, it dequeues the next invalidIndex from the 
	 * invalidIndicies queue and calls the handleBadValue function
	 * with it.
	 */
	void checkForOtherBadValues();

private:
	/**
	 * A Queue containing any invalid indicies that need to be
	 * delt with.
	 */
	QQueue<invalidIndex> invalidInicies; 
};



}
#endif //OPTIKA_TREEVIEW_HPP_
