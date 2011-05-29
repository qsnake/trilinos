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


#include "Sundance.hpp"
#include "SundanceDOFMapBuilder.hpp"
#include "SundanceRaviartThomas.hpp"
#include "SundanceFunctionSupportResolver.hpp"


int main(int argc, char** argv)
{
  int stat = 0;
  try
  {
    Sundance::init(&argc, &argv);

    /* */
    int npx = MPIComm::world().getNProc();

    /* Create a mesh. It will be of type BasisSimplicialMesh, and will
     * be built using a PartitionedLineMesher. */
    MeshType meshType = new BasicSimplicialMeshType();
    MeshSource mesher = new PartitionedRectangleMesher(0.0, 3.0, 2, npx,
      0.0, 1.0, 2, 1,
      meshType);
    Mesh mesh = mesher.getMesh();


    FieldWriter w = new VerboseFieldWriter();
    w.addMesh(mesh);
    w.write();

    FieldWriter w2 = new VTKWriter("rtMesh");
    w2.addMesh(mesh);
    w2.write();

    /* Create a cell filter that will identify the maximal cells
     * in the interior of the domain */
    CellFilter interior = new MaximalCellFilter();

    /* */
    Expr v = new TestFunction(new RaviartThomas(2));
    /* */
    Expr u = new UnknownFunction(new RaviartThomas(2));

    /* */
    Expr q = new TestFunction(new Lagrange(1));
    /* */
    Expr p = new UnknownFunction(new Lagrange(1));

    Out::os() << "v = " << describeFunction(v) << std::endl;
    Out::os() << "u = " << describeFunction(u) << std::endl;

    Out::os() << "q = " << describeFunction(q) << std::endl;
    Out::os() << "p = " << describeFunction(p) << std::endl;



    QuadratureFamily quad = new GaussianQuadrature(2);
    Expr eqn = Integral(interior, v*u + p*q, quad);
    Expr dum;

    RCP<FunctionSupportResolver> fsr 
      = rcp(new FunctionSupportResolver(eqn, dum, tuple(List(v,q).flatten()), tuple(List(u,p).flatten()),
          dum, dum, tuple(dum), false));
          
    
    int verb = 0;
    DOFMapBuilder builder(mesh, fsr, false, verb);

    for (int br = 0; br<builder.rowMap().size(); br++)
    {
      RCP<DOFMapBase> rm = builder.rowMap()[br];
      rm->print(Out::os());
    }

    
  }
	catch(std::exception& e)
  {
    stat = -1;
    std::cerr << "RT dof test FAILED" << std::endl;
    std::cerr << e.what() << std::endl;
  }

  return stat;
  
}
