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
#include "Optika_metawindow.hpp"
#include "Optika_Version.hpp"
#include <QLineEdit>
#include <QLabel>
#include <QPushButton>
#include <QFileDialog>
#include <QMessageBox> 
#include <QAction>
#include <QMenu>
#include <QMenuBar>
#include <QtGui>
#include <QIcon>
#include <iostream>
#include <algorithm>
namespace Optika{


const int numRecentDocuments = 7; 

SearchWidget::SearchWidget(TreeModel *treeModel, TreeView *treeView, QWidget *parent):
	QDialog(parent),
	treeModel(treeModel),
	treeView(treeView)
{
	matchesLabel = new QLabel(tr("Matches"));
	searchButton = new QPushButton(tr("Search"));
	connect(searchButton, SIGNAL(clicked(bool)), this, SLOT(search()));
	closeButton = new QPushButton(tr("Close"));
	connect(closeButton, SIGNAL(clicked(bool)), this, SLOT(hide()));
	nextButton = new QPushButton(tr("Next"));
	connect(nextButton, SIGNAL(clicked(bool)), this, SLOT(next()));
	previousButton = new QPushButton(tr("Previous"));
	connect(previousButton, SIGNAL(clicked(bool)), this, SLOT(previous()));
	searchTermsEdit = new QLineEdit(tr("Enter Search Terms Here"));
	QGridLayout *layout = new QGridLayout(this);
	layout->addWidget(searchTermsEdit,0,0);
	layout->addWidget(searchButton,0,1);
	layout->addWidget(nextButton,0,3);
	layout->addWidget(previousButton,0,2);
	layout->addWidget(closeButton,2,3);
	layout->addWidget(matchesLabel,3,0);
	setLayout(layout);
	nextButton->setDisabled(true);
	previousButton->setDisabled(true);
	setSizeGripEnabled(true);
	setWindowTitle(tr("Search..."));
}

void SearchWidget::search(){
	currentSearchResults = treeModel->match(treeModel->index(0,0,QModelIndex()), Qt::DisplayRole, searchTermsEdit->text(), 
	-1, Qt::MatchWrap | Qt::MatchContains | Qt::MatchRecursive);
	currentSearchResults = removeHiddenItems(currentSearchResults);
	currentSearchIterator = currentSearchResults.begin();	
	int searchSize = currentSearchResults.size();
	matchesLabel->setText("Matches ("+ QString::number(searchSize) + ")");
	if(searchSize <= 0){
		nextButton->setDisabled(true);
		previousButton->setDisabled(true);
	}
	else if(searchSize == 1){
		nextButton->setDisabled(true);
		previousButton->setDisabled(true);
		treeView->setCurrentIndex(*currentSearchIterator);
	}
	else{
		nextButton->setDisabled(false);
		previousButton->setDisabled(false);
		treeView->setCurrentIndex(*currentSearchIterator);
	}
}

void SearchWidget::next(){
	currentSearchIterator++;
	if(currentSearchIterator == currentSearchResults.end()){
		currentSearchIterator = currentSearchResults.begin();
	}
	treeView->setCurrentIndex(*currentSearchIterator);
}

void SearchWidget::previous(){
	currentSearchIterator--;
	if(currentSearchIterator == currentSearchResults.begin()-1){
		currentSearchIterator = currentSearchResults.end() -1;
	}
	treeView->setCurrentIndex(*currentSearchIterator);
}

QModelIndexList SearchWidget::removeHiddenItems(QModelIndexList& items){
	QModelIndexList toReturn;
	for(QModelIndexList::iterator it = items.begin(); it != items.end(); ++it){
		if(!treeView->isRowHidden(it->row(), it->parent())){
			toReturn.append(*it);
		}
	}
	return toReturn;
}

MetaWindow::MetaWindow(Teuchos::RCP<Teuchos::ParameterList> validParameters, QString fileName){
	model = new TreeModel(validParameters, fileName);
	initilization();
}

MetaWindow::MetaWindow(Teuchos::RCP<Teuchos::ParameterList> validParameters, void (*customFunc)(Teuchos::RCP<const Teuchos::ParameterList>), QString fileName){
	model = new TreeModel(validParameters, fileName);
	initilization(customFunc);
}

MetaWindow::MetaWindow(Teuchos::RCP<Teuchos::ParameterList> validParameters, Teuchos::RCP<DependencySheet> dependencySheet, QString fileName){
	model = new TreeModel(validParameters, dependencySheet, fileName);
	initilization();
} 

MetaWindow::MetaWindow(Teuchos::RCP<Teuchos::ParameterList> validParameters, Teuchos::RCP<DependencySheet> dependencySheet, void (*customFunc)(Teuchos::RCP<const Teuchos::ParameterList>), QString fileName){
	model = new TreeModel(validParameters, dependencySheet, fileName);
	initilization(customFunc);
} 

MetaWindow::~MetaWindow(){
	saveSettings();
}

void MetaWindow::closeEvent(QCloseEvent *event){
	if(!model->isSaved()){
		if(saveCurrentUnsavedFile()){
			saveSettings();
			event->accept();
		}
		else{
			event->ignore();
		}
	}
	else{
		saveSettings();
		event->accept();
	}
}

void MetaWindow::initilization(void (*customFunc)(Teuchos::RCP<const Teuchos::ParameterList>)){
	this->customFunc = customFunc;
	delegate = new Delegate;
	view = new TreeView(model, delegate);
	view->setEditTriggers(QAbstractItemView::DoubleClicked | QAbstractItemView::SelectedClicked);
	searchWidget = new SearchWidget(model, view, this);
	searchWidget->hide();
	QPushButton *submitButton = new QPushButton(tr("Submit"), this);
	connect(submitButton, SIGNAL(clicked(bool)), this, SLOT(submit()));
	QWidget *centerWidget = new QWidget(this);
	QGridLayout *centerWidgetLayout = new QGridLayout(centerWidget);
	centerWidgetLayout->addWidget(view,0,0);
	centerWidgetLayout->addWidget(submitButton,1,0,Qt::AlignRight);
	centerWidget->setLayout(centerWidgetLayout);
	setCentralWidget(centerWidget);
	createActions();
	createMenus();
	resize(800,600);
	currentLoadDir = QDir::homePath();
	currentSaveDir = QDir::homePath();
	loadLastSettings();
	setWindowTitle(tr("Parameter Input"));
	view->show();
	view->header()->resizeSections(QHeaderView::ResizeToContents);
	view->header()->setMovable(false);
}

void MetaWindow::createMenus(){
//	recentMenu = new QMenu(tr("Recent Solvers"));
//	QAction *noRecentAct = new QAction(tr("No Recent Documents"),this);
//	noRecentAct->setEnabled(false);
//	recentMenu->addAction(noRecentAct);
	fileMenu = menuBar()->addMenu(tr("File"));
	fileMenu->addAction(resetAct);
	//fileMenu->addMenu(recentMenu);
	fileMenu->addSeparator();
	fileMenu->addAction(saveAct);
	fileMenu->addAction(saveAsAct);
	fileMenu->addAction(loadAct);
	fileMenu->addSeparator();
	fileMenu->addAction(quitAct);
	helpMenu = menuBar()->addMenu(tr("Help"));
	helpMenu->addAction(aboutAct);
	helpMenu->addAction(searchAct);
}

void MetaWindow::createActions(){
	resetAct = new QAction(tr("&Reset"),this);
	resetAct->setShortcut(tr("Ctrl+R"));
	resetAct->setStatusTip(tr("Reset the list to its original state."));
	connect(resetAct, SIGNAL(triggered()), this, SLOT(resetModel()));

	saveAct = new QAction(tr("&Save"),this);
	saveAct->setShortcut(tr("Ctrl+S"));
	saveAct->setStatusTip(tr("Save the current file."));
	connect(saveAct, SIGNAL(triggered()), this, SLOT(saveFile()));

	saveAsAct = new QAction(tr("Save As..."),this);
	saveAsAct->setStatusTip(tr("Save the current file to a specified file name."));
	connect(saveAsAct, SIGNAL(triggered()), this, SLOT(saveFileAs()));

	loadAct = new QAction(tr("&Load"),this);
	loadAct->setShortcut(tr("Ctrl+L"));
	loadAct->setStatusTip(tr("Load input file"));
	connect(loadAct, SIGNAL(triggered()), this, SLOT(loadFile()));

	quitAct = new QAction(tr("&Quit"),this);
	quitAct->setShortcut(tr("Ctrl+Q"));
	quitAct->setStatusTip(tr("Quit"));
	connect(quitAct, SIGNAL(triggered()), this, SLOT(close()));

	aboutAct = new QAction(tr("About"),this);
	searchAct = new QAction(tr("Search For Parameter/Parameter List"), this);
	searchAct->setToolTip("Search for a particular Parameter or ParameterList");
	connect(aboutAct, SIGNAL(triggered()), this, SLOT(showAbout()));
	connect(searchAct, SIGNAL(triggered()), this, SLOT(initiateSearch()));
}

void MetaWindow::load(){
	QString fileName = QFileDialog::getOpenFileName(this, tr("Load..."), currentLoadDir, tr("Xml (*.xml)"));
	if(fileName != ""){
		model->readInput(fileName);
		currentLoadDir = fileName.section("/",0,-2);
		addRecentDocument(fileName);
	}
}

void MetaWindow::loadLastSettings(){
	QFile *file = new QFile("settings.xml");
	if(file->open(QIODevice::ReadOnly)){
		QXmlStreamReader xmlReader(file);
		while(!xmlReader.isEndDocument()){
			if(xmlReader.isStartElement()){
				if(xmlReader.name() == "lastsavedir"){
					QString tempCurSave = xmlReader.readElementText();
					if(tempCurSave != ""){
						currentSaveDir = tempCurSave;
					}
				}
				else if(xmlReader.name() == "lastloaddir"){
					QString tempCurLoad = xmlReader.readElementText();
					if(tempCurLoad != ""){
						currentLoadDir = tempCurLoad;
					}
				}
				else if(xmlReader.name() == "xres"){
					QString theWidth = xmlReader.readElementText();
					if(theWidth != ""){
						resize(theWidth.toInt(), height());
					}
				}
				else if(xmlReader.name() == "yres"){
					QString theHeight = xmlReader.readElementText();
					if(theHeight != ""){
						resize(width(), theHeight.toInt());
					}
				}
				else if(xmlReader.name() == "xpos"){
					QString xpos = xmlReader.readElementText();
					if(xpos != ""){
						move(xpos.toInt(), y());
					}
				}
				else if(xmlReader.name() == "ypos"){
					QString ypos = xmlReader.readElementText();
					if(ypos != ""){
						move(x(), ypos.toInt());
					}
				}
				else if(xmlReader.name() == "recentdoc"){
					addRecentDocument(xmlReader.readElementText());
				}
			}
			xmlReader.readNext();
		}
		file->close();
	}
	delete file;
}


void MetaWindow::saveSettings(){
	QFile *file = new QFile("settings.xml");
	file->open(QIODevice::WriteOnly);
	QXmlStreamWriter xmlWriter(file);

	xmlWriter.setAutoFormatting(true);
	xmlWriter.writeStartDocument();
	xmlWriter.writeStartElement("settings");
		xmlWriter.writeStartElement("lastsavedir");
		xmlWriter.writeCharacters(currentSaveDir);
		xmlWriter.writeEndElement();
		xmlWriter.writeStartElement("lastloaddir");
		xmlWriter.writeCharacters(currentLoadDir);
		xmlWriter.writeEndElement();
		xmlWriter.writeStartElement("xres");
		xmlWriter.writeCharacters(QString::number(width()));
		xmlWriter.writeEndElement();
		xmlWriter.writeStartElement("yres");
		xmlWriter.writeCharacters(QString::number(height()));
		xmlWriter.writeEndElement();
		xmlWriter.writeStartElement("xpos");
		xmlWriter.writeCharacters(QString::number(x()));
		xmlWriter.writeEndElement();
		xmlWriter.writeStartElement("ypos");
		xmlWriter.writeCharacters(QString::number(y()));
		xmlWriter.writeEndElement();
		for(int i =0; i<recentDocsList.size(); ++i){
			xmlWriter.writeStartElement("recentdoc");
				xmlWriter.writeCharacters(recentDocsList.at(i));
			xmlWriter.writeEndElement();
		}
	xmlWriter.writeEndElement();
	xmlWriter.writeEndDocument();

	file->close();
	delete file;
}
	

void MetaWindow::addRecentDocument(QString recentDocument){
	recentDocsList.prepend(recentDocument);
	if(recentDocsList.size() > numRecentDocuments){
		recentDocsList.removeLast();
	}
//	updateRecentDocsMenu();
}

void MetaWindow::updateRecentDocsMenu(){
	recentMenu->clear();
	for(int i=0; i<recentDocsList.size(); ++i){
		QAction *recentDocAct = new QAction(recentDocsList.at(i).section("/",-1,-1),this);
		connect(recentDocAct, SIGNAL(triggered()), this, SLOT(loadRecentDoc()));
		recentMenu->addAction(recentDocAct);
	}
}

void MetaWindow::resetModel(){
	if(!model->isSaved()){
		saveCurrentUnsavedFile();
	}
	model->reset();
}

bool MetaWindow::saveFileAs(){
	QString fileName = QFileDialog::getSaveFileName(this, tr("Save To..."), currentSaveDir, tr("XML (*.xml)"));
	if(fileName.toStdString() != ""){
		if(!fileName.endsWith(".xml")){
			fileName = fileName.append(".xml");
		}
		if(model->writeOutput(fileName)){
			currentSaveDir = fileName.section("/",0,-2);
			addRecentDocument(fileName);		
			return true;
		}
	}
	return false;
}

void MetaWindow::saveFile(){
	QString currentFileName = model->getSaveFileName();
	if(currentFileName != ""){
		model->writeOutput(currentFileName);
	}
	else{
		saveFileAs();
	}
}


void MetaWindow::loadFile(){
	if(!model->isSaved()){
		saveCurrentUnsavedFile();
	}
	load();
}

bool MetaWindow::saveCurrentUnsavedFile(){
		QMessageBox saveQuestion(QMessageBox::Question, 
		tr("Save?"),
		tr("These choices have not been saved since you last made changes. Would you like to save them now?"), 
		QMessageBox::Yes | QMessageBox::No,
		this);
		saveQuestion.setDefaultButton(QMessageBox::Yes);
		int shouldSave = saveQuestion.exec(); 
		if(shouldSave == QMessageBox::Yes){
			return saveFileAs();
		}
		return true;
}

void MetaWindow::loadRecentDoc(){
	QString docName = dynamic_cast<QAction*>(sender())->text();
	int i =0;
	for(; i<recentDocsList.size();++i){
		if(recentDocsList.at(i).contains(docName)){
			break;
		}
	}
	if(!model->isSaved()){
		if(saveCurrentUnsavedFile()){
			model->readInput(recentDocsList.at(i));
		}
	}
}

void MetaWindow::showAbout(){
	QString aboutString = aboutInfo + "\nThis input obtainer was generated by Kurtis Nusbaum's Optika package, part of the Trilinos Project.\n\nVersion: " + QString::fromStdString(Optika_Version()) + "\nWebsite: trilinos.sandia.gov/packages/optika\nLicense: LGPL\nContact: klnusbaum@gmail.com";
	QMessageBox::about(this,
	"Optika Input Obtainer\n",
	aboutString);
}

void MetaWindow::initiateSearch(){
	searchWidget->show();
}

void MetaWindow::submit(){
	if(customFunc == 0){
		close();
	}
	else{
		(*customFunc)(model->getCurrentParameters());	
	}
}

void MetaWindow::setAboutInfo(QString aboutInfo){
	this->aboutInfo = aboutInfo;
}

QString MetaWindow::getAboutInfo(){
	return aboutInfo;
}

}

