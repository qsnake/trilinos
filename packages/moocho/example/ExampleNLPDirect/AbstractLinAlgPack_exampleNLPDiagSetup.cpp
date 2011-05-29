// ////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_exampleNLPDiagSetup.hpp

#include <assert.h>

#include <fstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "AbstractLinAlgPack_exampleNLPDiagSetup.hpp"
#include "AbstractLinAlgPack_VectorSpaceSerial.hpp"
#include "OptionsFromStreamPack_OptionsFromStream.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#ifdef USE_EPETRA_THYRA

#include "TSFCoreSerialVectorSpaceStd.hpp"
#include "AbstractLinAlgPack_VectorSpaceTSFCore.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "Epetra_Map.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#endif // USE_EPETRA_THYRA

/** \brief . */
int AbstractLinAlgPack::exampleNLPDiagSetup(
  int argc, char* argv[], MPI_Comm comm
  ,Teuchos::RCP<const VectorSpace>   *vec_space
  ,int *n, value_type *xo, bool *has_bounds, bool *dep_bounded
  )
{

  using std::endl;
  using std::setw;
  namespace mmp = MemMngPack;
  using Teuchos::RCP;
  typedef AbstractLinAlgPack::size_type size_type;
  typedef AbstractLinAlgPack::value_type value_type;

  using AbstractLinAlgPack::VectorSpace;
  using AbstractLinAlgPack::Vector;
  using AbstractLinAlgPack::VectorMutable;

  using Teuchos::CommandLineProcessor;

  // Get an idea of what processors we have.
  int num_proc, proc_rank;
  MPI_Comm_size( comm, &num_proc );
  MPI_Comm_rank( comm, &proc_rank );

  // Get the size of the problem to solve
  *n = 4;
  // Get the starting point
  *xo = 0.1;
  // Determine if the NLP has bounds or not.
  *has_bounds = false;
  // Make the dependent or independent variables bounded.
  *dep_bounded = true;
#ifdef USE_EPETRA_THYRA
  // Serial or parallel?
  bool in_parallel = false;
  // Use TSF?
  bool use_tsf = false;
#endif // USE_EPETRA_THYRA

  CommandLineProcessor  command_line_processor;
  
  command_line_processor.setOption( "n",  n,   "Global number of dependent (and independent) variables" );
  command_line_processor.setOption( "xo", xo,  "Initial guess of the solution" );
  command_line_processor.setOption(
    "has-bounds", "no-has-bounds", has_bounds
    ,"Determine if the NLP has bounds or not" );
  command_line_processor.setOption(
    "dep-bounded", "indep-bounded", dep_bounded
    ,"Determine if the dependent or independent variables are bounded" );
#ifdef USE_EPETRA_THYRA
  command_line_processor.setOption(
    "in-parallel", "in-serial", &in_parallel
    ,"Determine if computations are performed in parallel or not" );
  command_line_processor.setOption(
    "use-tsf", "no-use-tsf", &use_tsf
    ,"Determine whether to use TSF vectors or not" );
#endif // USE_EPETRA_THYRA
  
  CommandLineProcessor::EParseCommandLineReturn
    parse_return = command_line_processor.parse(argc,argv,&std::cerr);
  
  if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
    return parse_return;

  // Create the vector space object to use.

#ifdef USE_EPETRA_THYRA

  using AbstractLinAlgPack::VectorSpaceTSFCore;

  if(in_parallel) {
    //
    // Use parallel vectors!
    //
    Teuchos::RCP<Epetra_Comm> comm;
#ifdef HAVE_MPI
    comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    comm = Teuchos::rcp(new Epetra_SerialComm());
#endif
    Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(*n,0,*comm));
    Teuchos::set_extra_data(comm, "comm", Teuchos::outArg(map));
    *vec_space = Teuchos::rcp(new VectorSpaceTSFCore(Teuchos::rcp(new TSFCore::EpetraVectorSpace(map))));
  }
  else {
    //
    // Use serial vectors
    //
    if( use_tsf ) {
      *vec_space = Teuchos::rcp(new VectorSpaceTSFCore(Teuchos::rcp(new TSFCore::SerialVectorSpaceStd<value_type>(*n))));
    }
    else {
      *vec_space = Teuchos::rcp(new AbstractLinAlgPack::VectorSpaceSerial(*n));
    }
  }

#else // USE_EPETRA_THYRA

  *vec_space = Teuchos::rcp(new AbstractLinAlgPack::VectorSpaceSerial(*n));

#endif // USE_EPETRA_THYRA
  
  return 0;
}
