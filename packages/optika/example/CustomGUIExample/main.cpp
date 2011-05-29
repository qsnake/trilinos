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
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Array.hpp"	
#include "Teuchos_Version.hpp"
#include "Optika_GUI.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerboseObject.hpp"

int main(int argc, char* argv[])
{
 /*
  * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              ATTENTION              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  * !!!!   PLEASE VIEW THE BASIC EXAMPLE FIRST BEFORE READING THIS EXAMPLE. IT PROVIDES FUNDAMENTAL    !!!! 
  * !!!!   KNOWLEDGE THAT WILL BE VERY HELPFUL IN UNDERSTANDING THIS EXAMPLE.                          !!!!
  * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  */ 

  /*
   * Defautls not good enought for you? That's understandable and why we've given
   * you a few ways to customize your GUI as you see fit. The purpose of this example is
   * to demonstrate those capabilities.
   */

  /* 
   * So let's setup a list of parameters to obtain. They could be anything really. We'll just reuse the
   * list we created in the Basic Example (which you've read already because you're smart and headed the gigantic
   * warning I put at the top of this example(.
   */
  Teuchos::RCP<Teuchos::ParameterList> My_List = Teuchos::RCP<Teuchos::ParameterList>(new Teuchos::ParameterList);

  My_List->set("Max Iters", 1550, "Determines the maximum number of iterations in the solver");
  My_List->set("Tolerance", 1e-10, "The tolerance used for the convergence check");
  
  Teuchos::RCP<Optika::StringValidator> solverValidator = 
     Teuchos::RCP<Optika::StringValidator>(new Optika::StringValidator(Teuchos::tuple<std::string>("GMRES", "CG", "TFQMR")));
  My_List->set( "Solver", "GMRES", "The type of solver to use.", solverValidator);

  Teuchos::Array<double> doubleArray( 10, 0.0 );
  My_List->set("Initial Guess", doubleArray, "The initial guess as a RCP to an array object.");

  Teuchos::ParameterList&
    Prec_List = My_List->sublist("Preconditioner",false,"Sublist that defines the preconditioner.");

  Prec_List.set("Type", "ILU", "The tpye of preconditioner to use");
  Prec_List.set("Drop Tolerance", 1e-3
                ,"The tolerance below which entries from the\n""factorization are left out of the factors.");

  /*
   * Now here's were things get a little differnent. Instead of just calling Optika::getInput(),
   * we're actually going to create and OptikaGUI object. It will be the vehical through which
   * we customize the GUI.
   */
   Optika::OptikaGUI myGUI(My_List);

  /*
   * Now we can start configuring what our GUI will look like. Let's start with the window
   * title.
   */
   myGUI.setWindowTitle("My Custom GUI");

  /*
   * We can set the information that will be displayed in the about dialog for the
   * gui.
   */
   myGUI.setAboutInfo("This is a little GUI I made.");

  /*
   * If you have an icon you'd like to use as the WindowIcon you can do that too.
   * just specify the path to the file containig the icon. Supported formats will
   * vary from system to system and your QT installation, but the following
   * should always work:
   	-BMP   -GIF  -JPG  -JPEG
	-MNG   -PNG  -PBM  -PGM
	-PPM   -TIFF -XBM  -XPM
	-SVG
   */
   myGUI.setWindowIcon("myIcon.png");

  /*
   * Now if you really wanna dive into your GUI design you can use
   * QT style sheets. Since optika is build on top of QT, you can
   * use QT style sheets to style the various widgets used by
   * Optika. The main widges Optika uses thay you'll probably 
   * wanna style are:
      -QTreeView -QDialog -QMainWindow
      -QPushButton -QSpinBox -QMenu
	  -QMenuBar
   * You might need to look at some of the Optika source code to really
   * get fine-grained styling control. Once your stylesheet is made,
   * all you need to do is specify the filepath to the Optika_GUI
   * object. Also note the style sheet provided in this example is 
   * exceptionally ugly.
   */
   myGUI.setStyleSheet("myStyleSheet.qss");

  /*
   * Now that we're all ready to go, we just call the exec funtion on
   * our Optika_GUI object.
   */
   myGUI.exec();

   /*
    * That's it! You can make even more awesome GUIs now!.
	*/
  return 0;
}

