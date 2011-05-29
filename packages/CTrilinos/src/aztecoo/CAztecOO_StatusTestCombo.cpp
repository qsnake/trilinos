
/*! @HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */


#include "CTrilinos_config.h"


#ifdef HAVE_CTRILINOS_AZTECOO


#include "CTrilinos_enums.h"
#include "CAztecOO_StatusTestCombo.h"
#include "CAztecOO_StatusTestCombo_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CAztecOO_StatusTest_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"


//
// Definitions from CAztecOO_StatusTestCombo.h
//


extern "C" {


CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestCombo_ID_t>(id);
}

CTrilinos_Universal_ID_t AztecOO_StatusTestCombo_Generalize ( 
  CT_AztecOO_StatusTestCombo_ID_t id )
{
    return CTrilinos::abstractType<CT_AztecOO_StatusTestCombo_ID_t>(id);
}

CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Create ( 
  CT_ComboType_E_t t )
{
    return CAztecOO::storeNewStatusTestCombo(new AztecOO_StatusTestCombo(
        (AztecOO_StatusTestCombo::ComboType) t));
}

CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Create_OneTest ( 
  CT_ComboType_E_t t, CT_AztecOO_StatusTest_ID_t aID )
{
    const Teuchos::RCP<AztecOO_StatusTest> a = CAztecOO::getStatusTest(aID);
    return CAztecOO::storeNewStatusTestCombo(new AztecOO_StatusTestCombo(
        (AztecOO_StatusTestCombo::ComboType) t, *a));
}

CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Create_TwoTests ( 
  CT_ComboType_E_t t, CT_AztecOO_StatusTest_ID_t aID, 
  CT_AztecOO_StatusTest_ID_t bID )
{
    const Teuchos::RCP<AztecOO_StatusTest> a = CAztecOO::getStatusTest(aID);
    const Teuchos::RCP<AztecOO_StatusTest> b = CAztecOO::getStatusTest(bID);
    return CAztecOO::storeNewStatusTestCombo(new AztecOO_StatusTestCombo(
        (AztecOO_StatusTestCombo::ComboType) t, *a, *b));
}

void AztecOO_StatusTestCombo_Destroy ( 
  CT_AztecOO_StatusTestCombo_ID_t * selfID )
{
    CAztecOO::removeStatusTestCombo(selfID);
}

CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_AddStatusTest ( 
  CT_AztecOO_StatusTestCombo_ID_t selfID, 
  CT_AztecOO_StatusTest_ID_t aID )
{
    const Teuchos::RCP<AztecOO_StatusTest> a = CAztecOO::getStatusTest(aID);
    return CAztecOO::storeStatusTestCombo(&( CAztecOO::getStatusTestCombo(
        selfID)->AddStatusTest(*a) ));
}

boolean AztecOO_StatusTestCombo_ResidualVectorRequired ( 
  CT_AztecOO_StatusTestCombo_ID_t selfID )
{
    return ((CAztecOO::getConstStatusTestCombo(
        selfID)->ResidualVectorRequired()) ? TRUE : FALSE);
}

CT_AztecOO_StatusType_E_t AztecOO_StatusTestCombo_CheckStatus ( 
  CT_AztecOO_StatusTestCombo_ID_t selfID, int CurrentIter, 
  CT_Epetra_MultiVector_ID_t CurrentResVectorID, 
  double CurrentResNormEst, boolean SolutionUpdated )
{
    const Teuchos::RCP<Epetra_MultiVector> CurrentResVector = 
        CEpetra::getMultiVector(CurrentResVectorID);
    return (CT_AztecOO_StatusType_E_t)( CAztecOO::getStatusTestCombo(
        selfID)->CheckStatus(CurrentIter, CurrentResVector.getRawPtr(), 
        CurrentResNormEst, ((SolutionUpdated) != FALSE ? true : false)) );
}

CT_AztecOO_StatusType_E_t AztecOO_StatusTestCombo_GetStatus ( 
  CT_AztecOO_StatusTestCombo_ID_t selfID )
{
    return (CT_AztecOO_StatusType_E_t)( CAztecOO::getConstStatusTestCombo(
        selfID)->GetStatus() );
}

CT_ComboType_E_t AztecOO_StatusTestCombo_GetComboType ( 
  CT_AztecOO_StatusTestCombo_ID_t selfID )
{
    return (CT_ComboType_E_t)( CAztecOO::getConstStatusTestCombo(
        selfID)->GetComboType() );
}


} // extern "C"




#endif /* HAVE_CTRILINOS_AZTECOO */


