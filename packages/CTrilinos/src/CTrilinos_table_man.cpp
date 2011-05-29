#include "CTrilinos_table_man.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_TableRepos.hpp"


extern "C" {


/*! Copies the RCP from one table into a second table. The new ID
 *  will be returned from the function. Both the old and the new
 *  IDs will need to be removed from the tables in order to destroy
 *  the object. */
CTrilinos_Universal_ID_t CT_Alias(CTrilinos_Universal_ID_t aid, CTrilinos_Table_ID_t new_table)
{
    return CTrilinos::tableRepos().alias(aid, new_table, true);
}

/*! Removes the RCP from one table and puts it in another. *aid will
 *  hold the new struct value afterward. Only the new RCP will need
 *  to be removed in order to destroy the object. */
void CT_Migrate(CTrilinos_Universal_ID_t *aid, CTrilinos_Table_ID_t new_table)
{
    CTrilinos_Universal_ID_t newid = CTrilinos::tableRepos().alias(*aid, new_table, false);
    *aid = newid;
}

/*! Checks to see if the underlying object referenced by a table
 *  entry is dynamic_cast'able to a given type (can be used to
 *  distinguish, e.g., an Epetra_SerialComm from an Epetra_MpiComm). */
boolean CT_TypeCheck(CTrilinos_Universal_ID_t aid, CTrilinos_Table_ID_t type)
{
    bool ret = CTrilinos::tableRepos().typeCheck(aid, type);
    return (ret ? TRUE : FALSE);
}


} // extern "C"
