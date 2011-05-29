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
#ifndef OPTIKA_GUI_HPP_
#define OPTIKA_GUI_HPP_
#include "Optika_DependencySheet.hpp"
#include "Optika_metawindow.hpp"
namespace Optika{
	/**
	 * Retreives the input for a Teuchos Parameter List using a GUI. Note the Parameter List will be edited.
	 * All user input will be stored in it.
	 *
	 * @param validParameters A list of parameters from which the users may specify values.
	 */
	void getInput(Teuchos::RCP<Teuchos::ParameterList> validParameters);

	/**
	 * Retreives the input for a Teuchos Parameter List using a GUI. Note the Parameter List will be edited.
	 * All user input will be stored in it. Also runs the function specified whenever the user hits the submit
	 * button.
	 *
	 * @param validParameters A list of parameters from which the users may specify values.
	 * @param customFunc Custom function to run whenever the user clicks the submit button.
	 */
	void getInput(Teuchos::RCP<Teuchos::ParameterList> validParameters, void (*customFunc)(Teuchos::RCP<const Teuchos::ParameterList>));

	/**
	 * Retreives the input for a Teuchos Parameter List using a GUI. Note the Parameter List will be edited.
	 * All user input will be stored in it.
	 *
	 * @param validParameters A list of parameters from which the users may specify values.
	 * @param dependencySheet A sheet listing any dependencies between parameters in the validParameters
	 * ParameterList.
	 */
	void getInput(Teuchos::RCP<Teuchos::ParameterList> validParameters, Teuchos::RCP<DependencySheet> dependencySheet);

	/**
	 * Retreives the input for a Teuchos Parameter List using a GUI. Note the Parameter List will be edited.
	 * All user input will be stored in it. Also runs the function specified whenever the user hits the submit
	 * button.
	 *
	 * @param validParameters A list of parameters from which the users may specify values.
	 * @param dependencySheet A sheet listing any dependencies between parameters in the validParameters
	 * ParameterList.
	 * @param customFunc Custom function to run whenever the user clicks the submit button.
	 */
	void getInput(Teuchos::RCP<Teuchos::ParameterList> validParameters, Teuchos::RCP<DependencySheet> dependencySheet, void (*customFunc)(Teuchos::RCP<const Teuchos::ParameterList>));

/**
 * A class that allows the user to create and customize their Optika GUI.
 */
class OptikaGUI{
public:
	/**
	 * Constructs an OptikaGUI object.
	 *
	 * @param validParameters A list of parameters from which the users may specify values.
	 */
	OptikaGUI(Teuchos::RCP<Teuchos::ParameterList> validParameters);

	/**
	 * Constructs an OptikaGUI object.
	 *
	 * @param validParameters A list of parameters from which the users may specify values.
	 * @param dependencySheet A sheet listing any dependencies between parameters in the validParameters
	 * ParameterList.
	 */
	OptikaGUI(Teuchos::RCP<Teuchos::ParameterList> validParameters, Teuchos::RCP<DependencySheet> dependencySheet);

	/**
	 * Runs the GUI and gets the user input.
	 */
	void exec();

	/**
	 * Adds the information specified to the about dialog of the GUI.
	 *
	 * @param aboutInfo Information to be added to the about dialog of the GUI.
	 */
	 void setAboutInfo(std::string aboutInfo);

	/**
	 * Sets the title of the GUI window that is displayed to the user.
	 *
	 * @param title A string containing what the title of the GUI window should be.
	 */
	void setWindowTitle(std::string title);

	/**
	 * Sets the window icon to the image specified in the filePath.
	 *
	 * @param filePath File path to the image that should be used as
	 *  the window icon.
	 */
	void setWindowIcon(std::string filePath);

	/**
	 * Sets the QT style sheet that should be used for the GUI.
	 *
	 * @param filePath File path to the QT style sheet to be used for
	 * the GUI.
	 */
	void setStyleSheet(std::string filePath);

	/**
	 * Sets the custom function to be used in the GUI. When ever the
	 * user hits submit, this function will be run.
	 *
	 * @param The custom function to be run whenever the user hits submit.
	 */
	void setCustomFunction(void (*customFunc)(Teuchos::RCP<const Teuchos::ParameterList>));

	/**
	 * Gets the window title.
	 * 
	 * @return A string containing the window title.
	 */
	std::string getWindowTitle();

	/**
	 * Gets the file path describing the location of the file
	 * being used for the window icon.
	 *
	 * @return The file path describing the location of the file
	 * being used for the window icon.
	 */
	std::string getWindowIcon();

	/**
	 * Gets the file path describing the location of the file
	 * being used as the QT Style Sheet.
	 *
	 * @return The file path describing the location of the file
	 * being used as the QT Style Sheet.
	 */
	std::string getStyleSheet();

	/**
	 * Gets the information to be added to the about dialog of the GUI.
	 *
	 * @return the information to be added to the about dialog of the GUI.
	 */
	 std::string getAboutInfo();

private:

	/**
	 * A list of parameters from which the users may specify values.
	 */
	Teuchos::RCP<Teuchos::ParameterList> validParameters;

	/**
	 * A sheet listing any dependencies between parameters in the validParameters
	 */
	Teuchos::RCP<DependencySheet> dependencySheet;

	/**
	 * A string containing the window title.
	 */
	std::string title;

	/**
	 * File path to the image that should be used as the window icon.
	 */
	std::string iconFilePath;

	/**
	 * File path to the QT style sheet to be used for the GUI.
	 */
	std::string styleSheetFilePath;

	/**
	 * Information to be added to the about dialog of the GUI.
	 */
	 std::string aboutInfo;

	/**
	 * The custom function to be run whenever the user hits submit.
	 */
	void (*customFunc)(Teuchos::RCP<const Teuchos::ParameterList>);
};

}

#endif //OPTIKA_GUI_HPP_
