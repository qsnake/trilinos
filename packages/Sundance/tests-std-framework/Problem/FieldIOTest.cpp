/* <Ignore> */
/***************************************************************************
 * Copyright (C) 2009
 * Kevin Long
 * Texas Tech University
 *
 * This library is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of the
 * License, or (at your option) any later version.
 *  
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *                                                                         
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA                     
 *  
 */
#include "Sundance.hpp"
#include "SundanceMeshIOUtils.hpp"

/* </Ignore> */


#if defined(HAVE_SUNDANCE_EXODUS) && defined(Trilinos_DATA_DIR)

int main(int argc, char** argv)
{
  try
  {
    Sundance::init(&argc, &argv);


    MPIComm world = MPIComm::world();
#ifdef HAVE_MPI
    MPIComm self = MPIComm(MPI_COMM_SELF);
#else 
    MPIComm self = world;
#endif

    std::string infile = "cyl-coarse";
    int numProc = world.getNProc();
    

#ifdef HAVE_SUNDANCE_CHACO
    Out::os() << "have chaco!" << std::endl;
    /* If in parallel, partition the mesh */
    if (world.getNProc() > 1)
    {
      /* Partition on the root only */
      if (world.getRank()==0)
      {
        Out::os() << "root processor is working on the partitioning" << std::endl;
        MeshType meshType = new BasicSimplicialMeshType();
        MeshSource mesher = new ExodusMeshReader(infile, meshType, self);

        RCP<SerialPartitionerBase> part 
          = rcp(new FileIOChacoPartitioner("part"));

        serialPartition(part, numProc, mesher, infile);
      }
      /* everyone else waits for the root processor to finish partitioning */
      MPIComm::world().synchronize();
    }

    /* Now do a readback test on the mesh */
    double err = readbackTester(infile, world);
#else
    double err = 0.0;
    if (world.getRank()==0)
    {
      Out::os() << "partitioning..." << std::endl;
      err = readbackTester(infile, self);
      Out::os() << "done partitioning..." << std::endl;
    }
    else
    {
      Out::os() << "waiting..." << std::endl;
    }
    Out::os() << "synching..." << std::endl;
    MPIComm::world().synchronize();
    Out::os() << "sharing erro..." << std::endl;
    MPIComm::world().bcast(&err, 1, MPIComm::INT, 0);
#endif
    
    /* Use a very tight tolerance because the readback should be
     * essentially exact */
    Sundance::passFailTest(err, 1.0e-15);
  }
	catch(std::exception& e)
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); return Sundance::testStatus(); 
}
    

#else // don't have exodus


int main(int argc, char** argv)
{
  Sundance::init(&argc, &argv);
  std::cout << "dummy FieldIOTest PASSED. Enable exodus to run the actual test" << std::endl;
  Sundance::finalize();
  return 0;
}


#endif

    
    

    
    
