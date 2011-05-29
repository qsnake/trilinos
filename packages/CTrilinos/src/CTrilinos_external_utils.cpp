#include "CTrilinos_external_utils.h"
#include "CTrilinos_test_utils.hpp"
#include "CTrilinos_exceptions.hpp"
#include "CEpetra_MpiComm.h"


#ifdef HAVE_MPI
#include "mpi.h"
#endif


extern "C" {


#ifdef HAVE_MPI

/*! Create an Epetra_MpiComm from Fortran */
CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Fortran_Create ( int fcomm )
{
    MPI_Fint mfcomm = (MPI_Fint) fcomm;
    MPI_Comm ccomm = MPI_Comm_f2c(mfcomm);

    /* duplicate the communicator so that we won't get any cross-talk
     * from the application */
    MPI_Comm dupcomm;
    int ret = MPI_Comm_dup(ccomm, &dupcomm);
    if (ret != MPI_SUCCESS)
        throw CTrilinos::CTrilinosMiscException("Error on MPI_Comm_dup");

    return Epetra_MpiComm_Create(dupcomm);
}

#endif /* HAVE_MPI */

/* Clear the tables between tests */
void CTrilinos_CleanSlate (  )
{
  CTrilinos::purgeAllTables();
}

} // extern "C"
