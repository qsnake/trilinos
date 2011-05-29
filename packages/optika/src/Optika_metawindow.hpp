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
#ifndef METAWINDOW_HPP_
#define METAWINDOW_HPP_
#include <QMainWindow>
#include <QDialog>
#include <QModelIndex>
#include "Optika_treeview.hpp"

class QAction;
class QMenu;
class QLabel;
class QPushButton;
class QLineEdit;
namespace Optika{

/**
 * A small widget that searchs through a parameter list for a particular name
 * of either a parameter or another parameter list.
 */
class SearchWidget : public QDialog{
	Q_OBJECT
public:
	/**
	 * Constructs a SearchWidget.
	 *
	 * @param treeModel The TreeModel being searched.
	 * @param treeView The TreeView being used to display the model.
	 * @param parent The parent widget.
	 */
	SearchWidget(TreeModel *treeModel, TreeView *treeView, QWidget *parent=0);

private slots:
	/**
	 * Searches the for a parameter or parameter list containing the string enterd
	 * in the search terms box.
	 */
	void search();

	/**
	 * Highlights the next result in the list of results that are set
	 * by the search function.
	 */
	void next();

	/**
	 * Highlights the previous result in the list of results that are set
	 * by the search function.
	 */
	void previous();

private:
	/**
	 * Removes any indicies in a QModelIndexList that correspond to a
	 * hidden item.
	 *
	 * @param items A list of indicies from which all hidden items
	 * will be removed
	 * @return A QModelIndexList identical to the one specified in the
	 * items parameter except that all indicies corresponding to hidden
	 * items have been removed.
	 */
	QModelIndexList removeHiddenItems(QModelIndexList& items);

	/**
	 * Widgets comprising a search widget
	 */
	QPushButton *searchButton, *closeButton, *nextButton, *previousButton;
	QLineEdit *searchTermsEdit;
	QLabel *matchesLabel;
	TreeModel *treeModel;
	TreeView *treeView;

	/**
	 * The results of the search last performed.
	 */
	QList<QModelIndex> currentSearchResults;

	/**
	 * An iterator over the results of the last search performed.
	 */
	QList<QModelIndex>::const_iterator currentSearchIterator;

};

/**
 * The Main Window that contains all other widgets in the Optika GUI.
 * For all undocumented functions please refer to the Qt API.
 */
class MetaWindow : public QMainWindow{
	Q_OBJECT
public:
	/**
	 * Constructs a MainWindow object.
	 * 
	 * @param validParameters The Parameter List the metawindow will display and the user will edit.
	 * @param fileName The name of a save file that may store previous values used by a user for the 
	 * Parameter List specified by validParameters.
	 */
	MetaWindow(Teuchos::RCP<Teuchos::ParameterList> validParameters, 
	QString fileName=QString());

	/**
	 * Constructs a MainWindow object.
	 * 
	 * @param validParameters The Parameter List the metawindow will display and the user will edit.
	 * @param customFunc The function to run whenever the user clicks the submit button.
	 * @param fileName The name of a save file that may store previous values used by a user for the 
	 * Parameter List specified by validParameters.
	 */
	MetaWindow(Teuchos::RCP<Teuchos::ParameterList> validParameters, 
	void (*customFunc)(Teuchos::RCP<const Teuchos::ParameterList>),
	QString fileName=QString());

	/**
	 * Constructs a MainWindow object.
	 * 
	 * @param validParameters The Parameter List the metawindow will display and the user will edit.
	 * @param dependencySheet A sheet listing any dependencies between parameters in the validParameters
	 * ParameterList.
	 * @param fileName The name of a save file that may store previous values used by a user for the 
	 * Parameter List specified by validParameters.
	 */
	MetaWindow(Teuchos::RCP<Teuchos::ParameterList> validParameters, 
	Teuchos::RCP<DependencySheet> dependencySheet,
	QString fileName=QString());

	/**
	 * Constructs a MainWindow object.
	 * 
	 * @param validParameters The Parameter List the metawindow will display and the user will edit.
	 * @param dependencySheet A sheet listing any dependencies between parameters in the validParameters
	 * ParameterList.
	 * @param customFunc The function to run whenever the user clicks the submit button.
	 * @param fileName The name of a save file that may store previous values used by a user for the 
	 * Parameter List specified by validParameters.
	 */
	MetaWindow(Teuchos::RCP<Teuchos::ParameterList> validParameters, 
	Teuchos::RCP<DependencySheet> dependencySheet,
	void (*customFunc)(Teuchos::RCP<const Teuchos::ParameterList>),
	QString fileName=QString());

	/**
	 * Deconstructer for the metawindow
	 */
	~MetaWindow();

	/**
	 * Adds the information specified to the about dialog of the GUI.
	 *
	 * @param aboutInfo Information to be added to the about dialog of the GUI.
	 */
	 void setAboutInfo(QString aboutInfo);

	/**
	 * Gets the information to be added to the about dialog of the GUI.
	 *
	 * @return the information to be added to the about dialog of the GUI.
	 */
	 QString getAboutInfo();


protected:
	/**
	 * Handles any QCloseEvents for the metawindow.
	 *
	 * @param event The QCloseEvent that was issued.
	 */
	void closeEvent(QCloseEvent *event);


private:
	/**
	 * Widgets comprising the MetaWindow
	 */
	SearchWidget *searchWidget;
	QAction *resetAct, *loadAct, *saveAct, *saveAsAct, *quitAct, *aboutAct, *searchAct;
	QMenu *fileMenu, *recentMenu, *helpMenu;

	/**
	 * Any additional about information that should be displayed in the about dialog.
	 */
	QString aboutInfo;

	/*
	 * The custom function to run when the user hits the submit button.
	 */
	void (*customFunc)(Teuchos::RCP<const Teuchos::ParameterList>);

	/**
	 * Load and save directory paths
	 */
	QString currentLoadDir, currentSaveDir;

	/**
	 * A list of recently used documents.
	 */
	QStringList recentDocsList;

	/**
	 * The TreeView being used in the metawindow.
	 */
	TreeView *view;

	/**
	 * The TreeModel being used to display the inputs.
	 */
	TreeModel *model;

	/**
	 * The deleages being used to modify any input values.
	 */
	Delegate *delegate;

	/**
	 * Common initialization shared by both constructors.
	 *
	 * @param customFunc The function to run whenever the user clicks the submit button.
	 */
	void initilization(void (*customFunc)(Teuchos::RCP<const Teuchos::ParameterList>)=0);

	/**
	 * Creates all the menus for the metawindow.
	 */
	void createMenus();

	/**
	 * Creates all necessary actions used in the menut items.
	 */
	void createActions();

	/**
	 * Loads previous parameter settings.
	 */
	void load();

	/**
	 * Loads the last state of the MetaWindow (things like window size and screen position).
	 */
	void loadLastSettings();

	/**
	 * Saves the state of the MetaWindow (things like window size and screen position).
	 */
	void saveSettings();

	/**
	 * Currently under developement
	 */
	void addRecentDocument(QString recentDocument);

	/**
	 * Currently under developement
	 */
	void updateRecentDocsMenu();

private slots:
	/**
	 * Resets the treemodel to its default state.
	 */
	void resetModel();

	/**
	 * Saves the parameter list settings to a user specified file.
	 */
	bool saveFileAs();

	/**
	 * Saves the current solver to the file the user has already specified.
	 */
	void saveFile();

	/**
	 * Loads a solver the user was previously working on and had saved.
	 */
	void loadFile();

	/**
	 * Asks the user whether or not they would like to currently save the file they are working on.
	 * Should be used when the user has modified the file and is about to perform an action that would cause those modifiation to be lost.
	 */
	bool saveCurrentUnsavedFile();

	/**
	 * Loads a document from the set of recent documents
	 */
	void loadRecentDoc();

	/**
	 * Shows information about the program.
	 */
	void showAbout();

	/**
	 * Starts a search for a parituclar Parameter or ParameterList.
	 */
	void initiateSearch();
	
	/**
	 * What should happen when the user clicks the submit button.
	 */
	void submit();
};



}
#endif /* METAWINDOW_HPP_ */
