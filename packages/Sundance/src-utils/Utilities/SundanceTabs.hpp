/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
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
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#ifndef SUNDANCE_TABS_H
#define SUNDANCE_TABS_H

#include "SundanceDefs.hpp"

namespace Sundance
{
/**
 * Tabbing utility for output. Constructing a new Tabs object automatically
 * increments the number of tabs to be written. When the Tabs object goes out
 * of scope, the original tabs level is restored. 
 *
 * The tab size and character can be specified through the setTabSize() and
 * setTabChar() methods, for example,
 * \code
 * Tabs::setTabChar('*');
 * Tabs::setTabSize(4);
 * \endcode
 * The tab character can be set on an object-by-object basis
 * through a constructor argument.
 *
 * By default, a header giving the depth of tabs is written to each line; this
 * can simplify scanning by eye for when a given tab level is reached. 
 * This header can be turned off by calling
 * \code
 * Tabs::showDepth() = false;
 * \endcode
 * 
 * Example: the code
 * \code
 * void f()
 * {
 *   Tabs tab;
 *   cout << tab << "in f()" << std::endl;
 *   g();
 *   cout << tab << "leaving f()" << std::endl;
 * }
 *
 * void g()
 * {
 *   Tabs tab0;
 *   cout << tab0 << "in g()" << std::endl;
 *   for (int i=0; i<3; i++)
 *     {
 *       Tabs tab1();
 *       cout << tab1 << "i=" << i << std::endl;
 *     }
 *   cout << tab0 << "leaving g()" << std::endl;
 * }
 * \endcode
 * writes the following output 
 * \code
 * [0]  in f()
 * [1]    in g()
 * [2]------i=0
 * [2]------i=1
 * [2]------i=2
 * [1]    leaving g()
 * [0]  leaving f()
 * \endcode 
 */
class Tabs
{
public:
  /** Constructor increments tab level */
  Tabs(bool jump=true);

  /** Destructor decrements tab level */
  ~Tabs();

  /** 
   * Print to stream. This method is usually not called directly, as
   * tabs will usually be written with the insertion operator
   */
  void print(std::ostream& os) const ;

  /** Change the tab size. Default is 2.  */
  static void setTabSize(int ts) {tabSize() = ts;}

  /** Indicate whether to print the tab depth as a header for each line. */
  static bool& showDepth() {static bool rtn = true; return rtn;}

private:
  /** */
  static int& tabLevel() {static int rtn = 0; return rtn;}

  /** */
  static int& tabSize() {static int rtn = 2; return rtn;}

  bool jump_;

  int myLevel_;
};
}

namespace Sundance
{
/** \relates Tabs stream insertion operator for tab */
inline std::ostream& operator<<(std::ostream& os, const Tabs& t)
{
  t.print(os);
  return os;
}
}

#endif
