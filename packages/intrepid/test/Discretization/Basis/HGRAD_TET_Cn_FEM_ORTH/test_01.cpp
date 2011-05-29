// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov) or
//                    Robert Kirby (robert.c.kirby@ttu.edu)
//
// ************************************************************************
// @HEADER


/** \file
    \brief  Unit test of orthogonal tetrahedral polynomial basis class
    \author Created by R. Kirby
*/

#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Intrepid_HGRAD_TET_Cn_FEM_ORTH.hpp"
#include "Intrepid_CubatureDirectTetDefault.hpp"
#include "Shards_CellTopology.hpp"
#include "Intrepid_PointTools.hpp"

#include <iostream>

using namespace Intrepid;

/** \brief Tests for orthogonal basis on tets.  Tests diagonality of mass matrices
    and does a code comparison to FIAT for values of derivatives
    \param argc [in] - number of command-line arguments
    \param argv [in] - command-line arguments
*/
int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);
  
  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
  
  *outStream \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|                           Unit Test OrthogonalBases                         |\n" \
    << "|                                                                             |\n" \
    << "|     1) Tests orthogonality of tetrahedral orthogonal basis                  |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
    << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n" \
    << "|                      Robert Kirby (robert.c.kirby@ttu.edu)                  |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";
  
  int errorFlag  = 0;

  const int deg = 3;

  Basis_HGRAD_TET_Cn_FEM_ORTH<double,FieldContainer<double> > myBasis( deg );
  const int polydim = myBasis.getCardinality();
  
  // First, get a reference quadrature rule

  CubatureDirectTetDefault<double,FieldContainer<double> > myCub(2*deg);
  FieldContainer<double> cubPts( myCub.getNumPoints() , 3 );
  FieldContainer<double> cubWts( myCub.getNumPoints() );

  myCub.getCubature( cubPts , cubWts );


  // Tabulate the basis functions at the cubature points
  FieldContainer<double> basisAtCubPts( polydim , myCub.getNumPoints() );

  myBasis.getValues( basisAtCubPts , cubPts , OPERATOR_VALUE );


  // Now let's compute the mass matrix
  for (int i=0;i<polydim;i++) {
    for (int j=i;j<polydim;j++) {
      double cur = 0.0;
      for (int k=0;k<myCub.getNumPoints();k++) {
        cur += cubWts(k) * basisAtCubPts( i , k ) * basisAtCubPts( j , k );
      }
      if (i != j && fabs( cur ) > 100. * INTREPID_TOL) {
        std::cout << "not diagonal" << i << " " << j << " " << fabs( cur ) << std::endl;
        errorFlag++;
      }
      if (i == j && fabs( cur ) <= 100. * INTREPID_TOL) {
        std::cout << "zero on diagonal" << i << " " << j << std::endl;
      }
    }
  }
  
  shards::CellTopology myTet_4( shards::getCellTopologyData< shards::Tetrahedron<4> >() );  
  const int np_lattice = PointTools::getLatticeSize( myTet_4 , deg , 0 );
  FieldContainer<double> lattice( np_lattice , 3 );
  PointTools::getLattice<double,FieldContainer<double> >( lattice , 
                                                          myTet_4 , 
                                                          deg , 
                                                          0 , 
                                                          POINTTYPE_EQUISPACED );        
                                
  FieldContainer<double> dBasisAtLattice( polydim , np_lattice , 3 );
  myBasis.getValues( dBasisAtLattice , lattice , OPERATOR_D1 );

  const double fiat_vals[] = {
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    5.477225575051661e+00, 2.738612787525831e+00, 2.738612787525831e+00,
    5.477225575051661e+00, 2.738612787525831e+00, 2.738612787525831e+00,
    5.477225575051661e+00, 2.738612787525831e+00, 2.738612787525831e+00,
    5.477225575051661e+00, 2.738612787525831e+00, 2.738612787525831e+00,
    5.477225575051661e+00, 2.738612787525831e+00, 2.738612787525831e+00,
    5.477225575051661e+00, 2.738612787525831e+00, 2.738612787525831e+00,
    5.477225575051661e+00, 2.738612787525831e+00, 2.738612787525831e+00,
    5.477225575051661e+00, 2.738612787525831e+00, 2.738612787525831e+00,
    5.477225575051661e+00, 2.738612787525831e+00, 2.738612787525831e+00,
    5.477225575051661e+00, 2.738612787525831e+00, 2.738612787525831e+00,
    5.477225575051661e+00, 2.738612787525831e+00, 2.738612787525831e+00,
    5.477225575051661e+00, 2.738612787525831e+00, 2.738612787525831e+00,
    5.477225575051661e+00, 2.738612787525831e+00, 2.738612787525831e+00,
    5.477225575051661e+00, 2.738612787525831e+00, 2.738612787525831e+00,
    5.477225575051661e+00, 2.738612787525831e+00, 2.738612787525831e+00,
    5.477225575051661e+00, 2.738612787525831e+00, 2.738612787525831e+00,
    5.477225575051661e+00, 2.738612787525831e+00, 2.738612787525831e+00,
    5.477225575051661e+00, 2.738612787525831e+00, 2.738612787525831e+00,
    5.477225575051661e+00, 2.738612787525831e+00, 2.738612787525831e+00,
    5.477225575051661e+00, 2.738612787525831e+00, 2.738612787525831e+00,
    0.000000000000000e+00, 4.743416490252569e+00, 1.581138830084190e+00,
    0.000000000000000e+00, 4.743416490252569e+00, 1.581138830084190e+00,
    0.000000000000000e+00, 4.743416490252568e+00, 1.581138830084190e+00,
    0.000000000000000e+00, 4.743416490252569e+00, 1.581138830084190e+00,
    0.000000000000000e+00, 4.743416490252569e+00, 1.581138830084190e+00,
    0.000000000000000e+00, 4.743416490252569e+00, 1.581138830084190e+00,
    0.000000000000000e+00, 4.743416490252569e+00, 1.581138830084190e+00,
    0.000000000000000e+00, 4.743416490252568e+00, 1.581138830084190e+00,
    0.000000000000000e+00, 4.743416490252568e+00, 1.581138830084190e+00,
    0.000000000000000e+00, 4.743416490252569e+00, 1.581138830084190e+00,
    0.000000000000000e+00, 4.743416490252569e+00, 1.581138830084190e+00,
    0.000000000000000e+00, 4.743416490252569e+00, 1.581138830084190e+00,
    0.000000000000000e+00, 4.743416490252569e+00, 1.581138830084190e+00,
    0.000000000000000e+00, 4.743416490252569e+00, 1.581138830084190e+00,
    0.000000000000000e+00, 4.743416490252569e+00, 1.581138830084190e+00,
    0.000000000000000e+00, 4.743416490252569e+00, 1.581138830084190e+00,
    0.000000000000000e+00, 4.743416490252568e+00, 1.581138830084190e+00,
    0.000000000000000e+00, 4.743416490252568e+00, 1.581138830084190e+00,
    0.000000000000000e+00, 4.743416490252568e+00, 1.581138830084190e+00,
    0.000000000000000e+00, 4.743416490252569e+00, 1.581138830084190e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 4.472135954999580e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 4.472135954999580e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 4.472135954999579e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 4.472135954999580e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 4.472135954999580e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 4.472135954999580e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 4.472135954999580e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 4.472135954999579e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 4.472135954999579e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 4.472135954999580e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 4.472135954999580e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 4.472135954999580e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 4.472135954999580e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 4.472135954999580e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 4.472135954999580e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 4.472135954999580e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 4.472135954999579e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 4.472135954999579e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 4.472135954999579e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 4.472135954999580e+00,
    -3.074085229787880e+01, -1.024695076595960e+01, -1.024695076595960e+01,
    -1.024695076595960e+01, 0.000000000000000e+00, 0.000000000000000e+00,
    1.024695076595960e+01, 1.024695076595960e+01, 1.024695076595960e+01,
    3.074085229787880e+01, 2.049390153191920e+01, 2.049390153191920e+01,
    -2.049390153191920e+01, -6.831300510639732e+00, -6.831300510639732e+00,
    -1.706460100885310e-15, 3.415650255319865e+00, 3.415650255319865e+00,
    2.049390153191919e+01, 1.366260102127946e+01, 1.366260102127946e+01,
    -1.024695076595960e+01, -3.415650255319866e+00, -3.415650255319866e+00,
    1.024695076595960e+01, 6.831300510639730e+00, 6.831300510639730e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    -2.049390153191920e+01, -6.831300510639732e+00, -6.831300510639732e+00,
    -1.706460100885310e-15, 3.415650255319865e+00, 3.415650255319865e+00,
    2.049390153191919e+01, 1.366260102127946e+01, 1.366260102127946e+01,
    -1.024695076595960e+01, -3.415650255319866e+00, -3.415650255319866e+00,
    1.024695076595960e+01, 6.831300510639730e+00, 6.831300510639730e+00,
    -1.706460100885310e-15, -5.688200336284365e-16, -5.688200336284365e-16,
    -1.024695076595960e+01, -3.415650255319866e+00, -3.415650255319866e+00,
    1.024695076595960e+01, 6.831300510639730e+00, 6.831300510639730e+00,
    -1.706460100885310e-15, -5.688200336284365e-16, -5.688200336284365e-16,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    -7.937253933193772e+00, -2.381176179958132e+01, -7.937253933193772e+00,
    -7.937253933193772e+00, -1.058300524425836e+01, -5.291502622129182e+00,
    -7.937253933193772e+00, 2.645751311064589e+00, -2.645751311064591e+00,
    -7.937253933193772e+00, 1.587450786638754e+01, 0.000000000000000e+00,
    5.291502622129181e+00, -1.058300524425837e+01, -4.406061034464155e-16,
    5.291502622129181e+00, 2.645751311064589e+00, 2.645751311064590e+00,
    5.291502622129181e+00, 1.587450786638754e+01, 5.291502622129181e+00,
    1.852025917745213e+01, 2.645751311064588e+00, 7.937253933193770e+00,
    1.852025917745213e+01, 1.587450786638754e+01, 1.058300524425836e+01,
    3.174901573277509e+01, 1.587450786638754e+01, 1.587450786638754e+01,
    -5.291502622129182e+00, -1.587450786638755e+01, -5.291502622129182e+00,
    -5.291502622129182e+00, -2.645751311064592e+00, -2.645751311064591e+00,
    -5.291502622129182e+00, 1.058300524425836e+01, -8.812122068928310e-16,
    7.937253933193772e+00, -2.645751311064592e+00, 2.645751311064591e+00,
    7.937253933193772e+00, 1.058300524425836e+01, 5.291502622129181e+00,
    2.116601048851673e+01, 1.058300524425836e+01, 1.058300524425836e+01,
    -2.645751311064591e+00, -7.937253933193774e+00, -2.645751311064591e+00,
    -2.645751311064591e+00, 5.291502622129181e+00, -4.406061034464155e-16,
    1.058300524425836e+01, 5.291502622129181e+00, 5.291502622129181e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    -6.480740698407860e+00, -3.240370349203930e+00, -2.268259244442751e+01,
    -6.480740698407860e+00, -3.240370349203930e+00, -9.721111047611791e+00,
    -6.480740698407860e+00, -3.240370349203930e+00, 3.240370349203929e+00,
    -6.480740698407860e+00, -3.240370349203930e+00, 1.620185174601965e+01,
    -6.480740698407860e+00, -3.240370349203930e+00, -1.620185174601965e+01,
    -6.480740698407860e+00, -3.240370349203930e+00, -3.240370349203932e+00,
    -6.480740698407860e+00, -3.240370349203930e+00, 9.721111047611787e+00,
    -6.480740698407860e+00, -3.240370349203930e+00, -9.721111047611791e+00,
    -6.480740698407860e+00, -3.240370349203930e+00, 3.240370349203929e+00,
    -6.480740698407860e+00, -3.240370349203930e+00, -3.240370349203930e+00,
    6.480740698407860e+00, 3.240370349203930e+00, -9.721111047611791e+00,
    6.480740698407860e+00, 3.240370349203930e+00, 3.240370349203929e+00,
    6.480740698407860e+00, 3.240370349203930e+00, 1.620185174601965e+01,
    6.480740698407860e+00, 3.240370349203930e+00, -3.240370349203930e+00,
    6.480740698407860e+00, 3.240370349203930e+00, 9.721111047611789e+00,
    6.480740698407860e+00, 3.240370349203930e+00, 3.240370349203929e+00,
    1.944222209522358e+01, 9.721111047611791e+00, 3.240370349203930e+00,
    1.944222209522358e+01, 9.721111047611791e+00, 1.620185174601965e+01,
    1.944222209522358e+01, 9.721111047611791e+00, 9.721111047611789e+00,
    3.240370349203930e+01, 1.620185174601965e+01, 1.620185174601965e+01,
    0.000000000000000e+00, -1.833030277982336e+01, -4.582575694955840e+00,
    0.000000000000000e+00, -1.833030277982336e+01, -4.582575694955840e+00,
    0.000000000000000e+00, -1.833030277982336e+01, -4.582575694955839e+00,
    0.000000000000000e+00, -1.833030277982336e+01, -4.582575694955840e+00,
    0.000000000000000e+00, -3.055050463303894e+00, 1.527525231651946e+00,
    0.000000000000000e+00, -3.055050463303894e+00, 1.527525231651946e+00,
    0.000000000000000e+00, -3.055050463303894e+00, 1.527525231651946e+00,
    0.000000000000000e+00, 1.222020185321557e+01, 7.637626158259732e+00,
    0.000000000000000e+00, 1.222020185321557e+01, 7.637626158259732e+00,
    0.000000000000000e+00, 2.749545416973504e+01, 1.374772708486752e+01,
    0.000000000000000e+00, -1.222020185321557e+01, -3.055050463303894e+00,
    0.000000000000000e+00, -1.222020185321557e+01, -3.055050463303894e+00,
    0.000000000000000e+00, -1.222020185321557e+01, -3.055050463303894e+00,
    0.000000000000000e+00, 3.055050463303893e+00, 3.055050463303893e+00,
    0.000000000000000e+00, 3.055050463303893e+00, 3.055050463303893e+00,
    0.000000000000000e+00, 1.833030277982336e+01, 9.165151389911680e+00,
    0.000000000000000e+00, -6.110100926607787e+00, -1.527525231651947e+00,
    0.000000000000000e+00, -6.110100926607787e+00, -1.527525231651947e+00,
    0.000000000000000e+00, 9.165151389911678e+00, 4.582575694955839e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, -5.612486080160912e+00, -1.309580085370879e+01,
    0.000000000000000e+00, -5.612486080160912e+00, -1.309580085370879e+01,
    0.000000000000000e+00, -5.612486080160911e+00, -1.309580085370879e+01,
    0.000000000000000e+00, -5.612486080160912e+00, -1.309580085370879e+01,
    0.000000000000000e+00, -5.612486080160912e+00, -1.870828693386971e+00,
    0.000000000000000e+00, -5.612486080160912e+00, -1.870828693386971e+00,
    0.000000000000000e+00, -5.612486080160912e+00, -1.870828693386971e+00,
    0.000000000000000e+00, -5.612486080160911e+00, 9.354143466934852e+00,
    0.000000000000000e+00, -5.612486080160911e+00, 9.354143466934852e+00,
    0.000000000000000e+00, -5.612486080160912e+00, 2.057911562725668e+01,
    0.000000000000000e+00, 5.612486080160912e+00, -5.612486080160912e+00,
    0.000000000000000e+00, 5.612486080160912e+00, -5.612486080160912e+00,
    0.000000000000000e+00, 5.612486080160912e+00, -5.612486080160912e+00,
    0.000000000000000e+00, 5.612486080160912e+00, 5.612486080160912e+00,
    0.000000000000000e+00, 5.612486080160912e+00, 5.612486080160912e+00,
    0.000000000000000e+00, 5.612486080160912e+00, 1.683745824048274e+01,
    0.000000000000000e+00, 1.683745824048273e+01, 1.870828693386970e+00,
    0.000000000000000e+00, 1.683745824048273e+01, 1.870828693386970e+00,
    0.000000000000000e+00, 1.683745824048273e+01, 1.309580085370879e+01,
    0.000000000000000e+00, 2.806243040080456e+01, 9.354143466934854e+00,
    0.000000000000000e+00, 0.000000000000000e+00, -1.322875655532295e+01,
    0.000000000000000e+00, 0.000000000000000e+00, -1.322875655532295e+01,
    0.000000000000000e+00, 0.000000000000000e+00, -1.322875655532295e+01,
    0.000000000000000e+00, 0.000000000000000e+00, -1.322875655532295e+01,
    0.000000000000000e+00, 0.000000000000000e+00, -1.322875655532295e+01,
    0.000000000000000e+00, 0.000000000000000e+00, -1.322875655532295e+01,
    0.000000000000000e+00, 0.000000000000000e+00, -1.322875655532295e+01,
    0.000000000000000e+00, 0.000000000000000e+00, -1.322875655532295e+01,
    0.000000000000000e+00, 0.000000000000000e+00, -1.322875655532295e+01,
    0.000000000000000e+00, 0.000000000000000e+00, -1.322875655532295e+01,
    0.000000000000000e+00, 0.000000000000000e+00, -8.812122068928310e-16,
    0.000000000000000e+00, 0.000000000000000e+00, -8.812122068928310e-16,
    0.000000000000000e+00, 0.000000000000000e+00, -8.812122068928310e-16,
    0.000000000000000e+00, 0.000000000000000e+00, -8.812122068928310e-16,
    0.000000000000000e+00, 0.000000000000000e+00, -8.812122068928310e-16,
    0.000000000000000e+00, 0.000000000000000e+00, -8.812122068928310e-16,
    0.000000000000000e+00, 0.000000000000000e+00, 1.322875655532295e+01,
    0.000000000000000e+00, 0.000000000000000e+00, 1.322875655532295e+01,
    0.000000000000000e+00, 0.000000000000000e+00, 1.322875655532295e+01,
    0.000000000000000e+00, 0.000000000000000e+00, 2.645751311064591e+01,
    9.524704719832526e+01, 2.381176179958132e+01, 2.381176179958132e+01,
    -1.058300524425836e+01, -1.322875655532295e+01, -1.322875655532295e+01,
    -1.058300524425837e+01, 2.645751311064586e+00, 2.645751311064586e+00,
    9.524704719832526e+01, 7.143528539874394e+01, 7.143528539874394e+01,
    4.233202097703347e+01, 1.058300524425837e+01, 1.058300524425837e+01,
    -1.058300524425836e+01, -5.291502622129184e+00, -5.291502622129184e+00,
    4.233202097703344e+01, 3.174901573277508e+01, 3.174901573277508e+01,
    1.058300524425837e+01, 2.645751311064592e+00, 2.645751311064592e+00,
    1.058300524425836e+01, 7.937253933193769e+00, 7.937253933193769e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    4.233202097703347e+01, 1.058300524425837e+01, 1.058300524425837e+01,
    -1.058300524425836e+01, -5.291502622129184e+00, -5.291502622129184e+00,
    4.233202097703344e+01, 3.174901573277508e+01, 3.174901573277508e+01,
    1.058300524425837e+01, 2.645751311064592e+00, 2.645751311064592e+00,
    1.058300524425836e+01, 7.937253933193769e+00, 7.937253933193769e+00,
    2.935026245019504e-31, 7.337565612548760e-32, 7.337565612548760e-32,
    1.058300524425837e+01, 2.645751311064592e+00, 2.645751311064592e+00,
    1.058300524425836e+01, 7.937253933193769e+00, 7.937253933193769e+00,
    2.935026245019504e-31, 7.337565612548760e-32, 7.337565612548760e-32,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    4.024922359499622e+01, 6.037383539249433e+01, 2.012461179749811e+01,
    1.341640786499874e+01, -1.565247584249853e+01, -2.236067977499790e+00,
    -1.341640786499874e+01, -2.906888370749726e+01, -1.565247584249853e+01,
    -4.024922359499622e+01, 2.012461179749811e+01, -2.012461179749811e+01,
    -3.577708763999664e+01, 8.944271909999161e+00, -8.944271909999159e+00,
    -2.979040983896728e-15, -4.472135954999583e+00, 4.472135954999577e+00,
    3.577708763999663e+01, 4.472135954999577e+01, 2.683281572999747e+01,
    -4.919349550499538e+01, -1.118033988749895e+01, -1.565247584249853e+01,
    4.919349550499536e+01, 3.801315561749642e+01, 3.354101966249684e+01,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    1.788854381999832e+01, 2.683281572999748e+01, 8.944271909999161e+00,
    1.489520491948364e-15, -1.341640786499874e+01, -4.472135954999580e+00,
    -1.788854381999832e+01, 8.944271909999145e+00, -8.944271909999161e+00,
    -2.236067977499789e+01, -2.236067977499788e+00, -6.708203932499368e+00,
    2.236067977499789e+01, 2.012461179749810e+01, 1.565247584249852e+01,
    -8.937122951690183e-15, -2.979040983896728e-15, -2.979040983896728e-15,
    4.472135954999580e+00, 6.708203932499371e+00, 2.236067977499790e+00,
    -4.472135954999580e+00, 2.236067977499787e+00, -2.236067977499790e+00,
    -4.468561475845092e-15, -1.489520491948364e-15, -1.489520491948364e-15,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    3.485685011586675e+01, 1.161895003862225e+01, 5.809475019311125e+01,
    1.161895003862225e+01, 0.000000000000000e+00, -1.549193338482967e+01,
    -1.161895003862225e+01, -1.161895003862225e+01, -2.711088342345191e+01,
    -3.485685011586675e+01, -2.323790007724450e+01, 2.323790007724450e+01,
    2.323790007724450e+01, 7.745966692414833e+00, 2.840187787218773e+01,
    1.934943878227166e-15, -3.872983346207416e+00, -1.420093893609386e+01,
    -2.323790007724450e+01, -1.549193338482966e+01, 5.163977794943212e+00,
    1.161895003862225e+01, 3.872983346207417e+00, 9.036961141150641e+00,
    -1.161895003862225e+01, -7.745966692414831e+00, -2.581988897471612e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    -3.872983346207416e+01, -1.290994448735805e+01, 7.745966692414837e+00,
    -3.224906463711944e-15, 6.454972243679026e+00, -3.872983346207421e+00,
    3.872983346207415e+01, 2.581988897471610e+01, 4.647580015448898e+01,
    -1.936491673103708e+01, -6.454972243679027e+00, -1.290994448735804e+00,
    1.936491673103708e+01, 1.290994448735805e+01, 1.807392228230127e+01,
    -3.224906463711944e-15, -1.074968821237314e-15, -1.074968821237314e-15,
    -5.034878350069641e+01, -1.678292783356547e+01, -1.161895003862225e+01,
    5.034878350069641e+01, 3.356585566713093e+01, 3.872983346207415e+01,
    -8.384756805651052e-15, -2.794918935217017e-15, -2.794918935217017e-15,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    1.039230484541326e+01, 6.754998149518620e+01, 1.558845726811989e+01,
    1.039230484541326e+01, 2.598076211353317e+01, 8.660254037844389e+00,
    1.039230484541326e+01, -1.558845726811989e+01, 1.732050807568879e+00,
    1.039230484541326e+01, -5.715767664977296e+01, -5.196152422706632e+00,
    -6.928203230275511e+00, -1.039230484541326e+01, -1.039230484541326e+01,
    -6.928203230275511e+00, -3.464101615137756e+00, -3.464101615137756e+00,
    -6.928203230275511e+00, 3.464101615137752e+00, 3.464101615137753e+00,
    2.424871130596427e+01, -1.558845726811990e+01, 1.732050807568872e+00,
    2.424871130596427e+01, 3.983716857408416e+01, 2.251666049839540e+01,
    1.039230484541326e+02, 5.196152422706632e+01, 5.196152422706632e+01,
    4.618802153517007e+00, 3.002221399786055e+01, 6.928203230275511e+00,
    4.618802153517007e+00, 2.309401076758505e+00, 2.309401076758504e+00,
    4.618802153517007e+00, -2.540341184434353e+01, -2.309401076758502e+00,
    1.154700538379250e+00, -9.814954576223640e+00, -4.041451884327382e+00,
    1.154700538379250e+00, 1.096965511460288e+01, 5.196152422706630e+00,
    4.618802153517006e+01, 2.309401076758502e+01, 2.309401076758503e+01,
    1.154700538379252e+00, 7.505553499465138e+00, 1.732050807568878e+00,
    1.154700538379252e+00, -6.350852961085883e+00, -5.773502691896256e-01,
    1.154700538379251e+01, 5.773502691896255e+00, 5.773502691896256e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    9.000000000000000e+00, 2.700000000000000e+01, 4.500000000000000e+01,
    9.000000000000000e+00, 1.200000000000000e+01, 1.800000000000000e+01,
    9.000000000000000e+00, -2.999999999999998e+00, -8.999999999999996e+00,
    9.000000000000000e+00, -1.800000000000000e+01, -3.600000000000000e+01,
    -6.000000000000000e+00, 1.200000000000000e+01, -1.600000000000000e+01,
    -6.000000000000000e+00, -2.999999999999998e+00, -3.000000000000001e+00,
    -6.000000000000000e+00, -1.800000000000000e+01, 9.999999999999996e+00,
    -2.100000000000000e+01, -2.999999999999997e+00, -3.700000000000000e+01,
    -2.100000000000000e+01, -1.800000000000000e+01, 1.599999999999999e+01,
    -3.600000000000000e+01, -1.800000000000000e+01, -1.800000000000000e+01,
    -1.000000000000000e+01, -3.000000000000001e+01, 6.000000000000003e+00,
    -1.000000000000000e+01, -5.000000000000002e+00, -4.999999999999999e+00,
    -1.000000000000000e+01, 1.999999999999999e+01, -1.600000000000000e+01,
    1.500000000000000e+01, -5.000000000000002e+00, -7.000000000000002e+00,
    1.500000000000000e+01, 1.999999999999999e+01, 2.199999999999999e+01,
    3.999999999999999e+01, 1.999999999999999e+01, 1.999999999999999e+01,
    -1.300000000000000e+01, -3.900000000000001e+01, -8.999999999999998e+00,
    -1.300000000000000e+01, 2.599999999999999e+01, -4.000000000000002e+00,
    5.199999999999999e+01, 2.599999999999999e+01, 2.599999999999999e+01,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    7.348469228349536e+00, 3.674234614174768e+00, 5.511351921262152e+01,
    7.348469228349536e+00, 3.674234614174768e+00, 2.082066281365702e+01,
    7.348469228349536e+00, 3.674234614174768e+00, -1.347219358530748e+01,
    7.348469228349536e+00, 3.674234614174768e+00, -4.776504998427198e+01,
    7.348469228349536e+00, 3.674234614174768e+00, 3.796709101313927e+01,
    7.348469228349536e+00, 3.674234614174768e+00, 3.674234614174771e+00,
    7.348469228349536e+00, 3.674234614174768e+00, -3.061862178478972e+01,
    7.348469228349536e+00, 3.674234614174768e+00, 2.082066281365702e+01,
    7.348469228349536e+00, 3.674234614174768e+00, -1.347219358530748e+01,
    7.348469228349536e+00, 3.674234614174768e+00, 3.674234614174768e+00,
    -4.082482904638631e+00, -2.041241452319316e+00, -1.347219358530748e+01,
    -4.082482904638631e+00, -2.041241452319316e+00, -2.041241452319317e+00,
    -4.082482904638631e+00, -2.041241452319316e+00, 9.389710680668845e+00,
    -4.082482904638631e+00, -2.041241452319316e+00, -7.756717518813399e+00,
    -4.082482904638631e+00, -2.041241452319316e+00, 3.674234614174765e+00,
    -4.082482904638631e+00, -2.041241452319316e+00, -2.041241452319317e+00,
    3.021037349432586e+01, 1.510518674716293e+01, -1.347219358530748e+01,
    3.021037349432586e+01, 1.510518674716293e+01, 4.368256707963334e+01,
    3.021037349432586e+01, 1.510518674716293e+01, 1.510518674716292e+01,
    1.102270384252430e+02, 5.511351921262151e+01, 5.511351921262151e+01,
    0.000000000000000e+00, 4.500000000000001e+01, 9.000000000000000e+00,
    0.000000000000000e+00, 4.500000000000001e+01, 9.000000000000000e+00,
    0.000000000000000e+00, 4.499999999999999e+01, 9.000000000000000e+00,
    0.000000000000000e+00, 4.500000000000001e+01, 9.000000000000000e+00,
    0.000000000000000e+00, -1.000000000000000e+01, -5.999999999999999e+00,
    0.000000000000000e+00, -1.000000000000000e+01, -5.999999999999999e+00,
    0.000000000000000e+00, -1.000000000000000e+01, -5.999999999999999e+00,
    0.000000000000000e+00, 4.999999999999996e+00, 8.999999999999996e+00,
    0.000000000000000e+00, 4.999999999999996e+00, 8.999999999999996e+00,
    0.000000000000000e+00, 9.000000000000000e+01, 5.400000000000000e+01,
    0.000000000000000e+00, 2.000000000000000e+01, 4.000000000000000e+00,
    0.000000000000000e+00, 2.000000000000000e+01, 4.000000000000000e+00,
    0.000000000000000e+00, 2.000000000000000e+01, 4.000000000000000e+00,
    0.000000000000000e+00, -5.000000000000002e+00, -1.000000000000001e+00,
    0.000000000000000e+00, -5.000000000000002e+00, -1.000000000000001e+00,
    0.000000000000000e+00, 4.000000000000000e+01, 2.400000000000000e+01,
    0.000000000000000e+00, 5.000000000000000e+00, 1.000000000000000e+00,
    0.000000000000000e+00, 5.000000000000000e+00, 1.000000000000000e+00,
    0.000000000000000e+00, 9.999999999999996e+00, 5.999999999999998e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, 2.078460969082653e+01, 2.598076211353316e+01,
    0.000000000000000e+00, 2.078460969082653e+01, 2.598076211353316e+01,
    0.000000000000000e+00, 2.078460969082652e+01, 2.598076211353316e+01,
    0.000000000000000e+00, 2.078460969082653e+01, 2.598076211353316e+01,
    0.000000000000000e+00, 3.464101615137756e+00, -1.327905619136139e+01,
    0.000000000000000e+00, 3.464101615137756e+00, -1.327905619136139e+01,
    0.000000000000000e+00, 3.464101615137756e+00, -1.327905619136139e+01,
    0.000000000000000e+00, -1.385640646055102e+01, -6.350852961085886e+00,
    0.000000000000000e+00, -1.385640646055102e+01, -6.350852961085886e+00,
    0.000000000000000e+00, -3.117691453623979e+01, 4.676537180435969e+01,
    0.000000000000000e+00, -2.309401076758503e+01, 3.464101615137756e+00,
    0.000000000000000e+00, -2.309401076758503e+01, 3.464101615137756e+00,
    0.000000000000000e+00, -2.309401076758503e+01, 3.464101615137756e+00,
    0.000000000000000e+00, 5.773502691896256e+00, 1.154700538379249e+00,
    0.000000000000000e+00, 5.773502691896256e+00, 1.154700538379249e+00,
    0.000000000000000e+00, 3.464101615137754e+01, 4.503332099679081e+01,
    0.000000000000000e+00, -3.002221399786054e+01, -5.196152422706631e+00,
    0.000000000000000e+00, -3.002221399786054e+01, -5.196152422706631e+00,
    0.000000000000000e+00, 4.503332099679080e+01, 2.944486372867091e+01,
    0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
    0.000000000000000e+00, 6.363961030678929e+00, 3.181980515339464e+01,
    0.000000000000000e+00, 6.363961030678929e+00, 3.181980515339464e+01,
    0.000000000000000e+00, 6.363961030678928e+00, 3.181980515339464e+01,
    0.000000000000000e+00, 6.363961030678929e+00, 3.181980515339464e+01,
    0.000000000000000e+00, 6.363961030678929e+00, 2.121320343559643e+00,
    0.000000000000000e+00, 6.363961030678929e+00, 2.121320343559643e+00,
    0.000000000000000e+00, 6.363961030678929e+00, 2.121320343559643e+00,
    0.000000000000000e+00, 6.363961030678928e+00, -2.757716446627535e+01,
    0.000000000000000e+00, 6.363961030678928e+00, -2.757716446627535e+01,
    0.000000000000000e+00, 6.363961030678929e+00, -5.727564927611036e+01,
    0.000000000000000e+00, -3.535533905932738e+00, -7.778174593052023e+00,
    0.000000000000000e+00, -3.535533905932738e+00, -7.778174593052023e+00,
    0.000000000000000e+00, -3.535533905932738e+00, -7.778174593052023e+00,
    0.000000000000000e+00, -3.535533905932738e+00, 2.121320343559642e+00,
    0.000000000000000e+00, -3.535533905932738e+00, 2.121320343559642e+00,
    0.000000000000000e+00, -3.535533905932738e+00, 1.202081528017131e+01,
    0.000000000000000e+00, 2.616295090390225e+01, -7.778174593052024e+00,
    0.000000000000000e+00, 2.616295090390225e+01, -7.778174593052024e+00,
    0.000000000000000e+00, 2.616295090390225e+01, 4.171930009000628e+01,
    0.000000000000000e+00, 9.545941546018391e+01, 3.181980515339464e+01,
    0.000000000000000e+00, 0.000000000000000e+00, 2.700000000000001e+01,
    0.000000000000000e+00, 0.000000000000000e+00, 2.700000000000001e+01,
    0.000000000000000e+00, 0.000000000000000e+00, 2.700000000000000e+01,
    0.000000000000000e+00, 0.000000000000000e+00, 2.700000000000001e+01,
    0.000000000000000e+00, 0.000000000000000e+00, 2.700000000000001e+01,
    0.000000000000000e+00, 0.000000000000000e+00, 2.700000000000001e+01,
    0.000000000000000e+00, 0.000000000000000e+00, 2.700000000000001e+01,
    0.000000000000000e+00, 0.000000000000000e+00, 2.700000000000000e+01,
    0.000000000000000e+00, 0.000000000000000e+00, 2.700000000000000e+01,
    0.000000000000000e+00, 0.000000000000000e+00, 2.700000000000001e+01,
    0.000000000000000e+00, 0.000000000000000e+00, -8.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, -8.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, -8.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, -8.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, -8.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, -8.000000000000000e+00,
    0.000000000000000e+00, 0.000000000000000e+00, 1.299999999999999e+01,
    0.000000000000000e+00, 0.000000000000000e+00, 1.299999999999999e+01,
    0.000000000000000e+00, 0.000000000000000e+00, 1.299999999999999e+01,
    0.000000000000000e+00, 0.000000000000000e+00, 9.000000000000000e+01
  };

  int fiat_index_cur = 0;
  for (int i=0;i<polydim;i++) {
    for (int j=0;j<np_lattice;j++) {
      for (int k=0;k<3;k++) {
        if (std::abs( dBasisAtLattice(i,j,k) - fiat_vals[fiat_index_cur] ) > 10.0*INTREPID_TOL ) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          
          // Output the multi-index of the value where the error is:
          *outStream << " At multi-index { ";
          *outStream << i << " " << j << " " << k;
          *outStream << "}  computed value: " << dBasisAtLattice(i,j,k)
                     << " but correct value: " << fiat_vals[fiat_index_cur] << "\n";
          *outStream << "Difference: " << std::abs( dBasisAtLattice(i,j,k) - fiat_vals[fiat_index_cur] ) << "\n";
        }
        fiat_index_cur++;
      }
    }
  }


  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";
  
  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
  
  return errorFlag;
}