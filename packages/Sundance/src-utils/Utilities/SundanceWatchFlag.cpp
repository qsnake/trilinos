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

#include "SundanceWatchFlag.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;


WatchFlag::WatchFlag(const std::string& name,
  const ParameterList& params)
  : name_(name), params_(rcp(new ParameterList(params)))
{
  if (name_.size() > 0) isActiveMap().put(name, true);
  else isActiveMap().put(name, false);
}

void WatchFlag::activate() 
{
  isActiveMap()[name()] = true;
}

void WatchFlag::deactivate() 
{
  isActiveMap()[name()] = false;
}

bool WatchFlag::isActive() const 
{
  return isActiveMap().get(name());
}

XMLObject WatchFlag::toXML() const 
{
  XMLObject xml("WatchFlag");
  xml.addAttribute("name", name());
  return xml;
}

int WatchFlag::param(const std::string& name) const 
{
  return params_->get<int>(name);
}


void WatchFlag::setParam(const std::string& name, int val)
{
  params_->set<int>(name, val);
}




RCP<ParameterList> WatchFlag::defaultParams()
{
  static RCP<ParameterList> rtn=rcp(new ParameterList());
  static bool first=true;
  if (first)
  {
    rtn->set<int>("evaluation", 0);
    rtn->set<int>("discrete function evaluation", 0);
    rtn->set<int>("symbolic preprocessing", 0);
    rtn->set<int>("equation set setup", 0);
    rtn->set<int>("assembler setup", 0);
    rtn->set<int>("assembly loop", 0);
    rtn->set<int>("dof map setup", 0);
    rtn->set<int>("dof map access", 0);
    rtn->set<int>("eval mediator", 0);
    rtn->set<int>("integration setup", 0);
    rtn->set<int>("integration", 0);
    rtn->set<int>("integral transformation", 0);    
    rtn->set<int>("fill", 0);
    rtn->set<int>("matrix config", 0);
    rtn->set<int>("vector config", 0);
    rtn->set<int>("solve control", 0);
    first = false;
  }
  return rtn;
}
