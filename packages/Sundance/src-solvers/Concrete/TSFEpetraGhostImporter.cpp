/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// **********************************************************************/

#include "TSFEpetraGhostImporter.hpp"
#include "TSFEpetraGhostView.hpp"
#include "TSFEpetraVector.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorImpl.hpp"
#endif

using namespace Teuchos;
using namespace TSFExtended;

EpetraGhostImporter
::EpetraGhostImporter(const RCP<const Epetra_Map>& localMap,
                      int nGhost,
                      const int* ghostElements)
  : localMap_(localMap),
    ghostMap_(),
    importer_()
{
  if (false && nGhost==0)
    {
      ghostMap_ = localMap_;

      importer_ = rcp(new Epetra_Import(*ghostMap_, *localMap_));
    }
  else
    {
      //bvbw not used      int nGlobal = localMap_->NumGlobalElements();
      int nLocal = localMap_->NumMyElements();
      int nGhostView = nLocal+nGhost;
      std::vector<int> globalIndices(nGhostView);
      for (int i=0; i<nLocal; i++) globalIndices[i] = localMap_->GID(i);
      for (int i=0; i<nGhost; i++) globalIndices[i+nLocal] = ghostElements[i];

      const Epetra_Comm& comm = localMap_->Comm();

      ghostMap_ = rcp(new Epetra_Map(-1, nGhostView, 
                                     &(globalIndices[0]), 0, comm));

      importer_ = rcp(new Epetra_Import(*ghostMap_, *localMap_));
    }
}

void EpetraGhostImporter
::importView(const Vector<double>& x,
             RCP<GhostView<double> >& ghostView) const
{
  /* If given an uninitialized ghost view, create a EpetraGhostView */
  if (ghostView.get()==0) 
    {
      ghostView = rcp(new EpetraGhostView());
    }

  /* Ensure that the ghost view contains an EpetraGhostView */
  EpetraGhostView* epgv 
    = dynamic_cast<EpetraGhostView*>(ghostView.get());

  TEST_FOR_EXCEPTION(epgv==0, std::runtime_error,
                     "argument ghostView to EpetraGhostImporter::importView() "
                     "could not be cast to a EpetraGhostView pointer");

  const Epetra_Vector& xVec = EpetraVector::getConcrete(x);

  /* Do the import */
  epgv->import(*importer_, xVec);
}
    
