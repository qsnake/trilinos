// @HEADER
// ************************************************************************
// 
//        TriKota: A Trilinos Wrapper for the Dakota Framework
//                  Copyright (2009) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// 
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include <iostream>
#include "TriKota_Driver.hpp"
#include "Teuchos_VerboseObject.hpp"

using namespace std;
using namespace Dakota;

// Dakota driver when linking in as a library 
// Assumes MPI_COMM_WORLD both for Dakota and the model evaluation
TriKota::Driver::Driver(const char* dakota_in,  
                               const char* dakota_out,
                               const char* dakota_err,
                               const char* dakota_restart_out)
 :  parallel_lib(), problem_db(parallel_lib)
{

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream(); 

  *out << "\nStarting TriKota_Driver!" << endl;

  parallel_lib.specify_outputs_restart(dakota_out, dakota_err, NULL,
                                       dakota_restart_out, 0);
  problem_db.manage_inputs(dakota_in);

  // instantiate the strategy
  selected_strategy = Strategy(problem_db);
}

MPI_Comm TriKota::Driver::getAnalysisComm()
{
  Model& first_model = *(problem_db.model_list().begin());
  MPI_Comm analysis_comm =
     first_model.parallel_configuration_iterator()->ea_parallel_level().server_intra_communicator();

  return analysis_comm;
}

ProblemDescDB& TriKota::Driver::getProblemDescDB()
{
  return problem_db;
}
  
void TriKota::Driver::run(Dakota::DirectApplicInterface* appInterface)
{

  Model& first_model = *(problem_db.model_list().begin());
  Interface& interface  = first_model.interface();

  // Pass a pointer to a Dakota::DirectApplicInterface
  interface.assign_rep(appInterface, false);

  selected_strategy.run_strategy();
}

const Dakota::Variables TriKota::Driver::getFinalSolution() const
{
  return selected_strategy.variables_results();
}
