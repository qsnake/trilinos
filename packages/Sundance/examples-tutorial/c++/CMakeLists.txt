

INCLUDE(AddTestBatch)

SET(MPITests
  PoissonDemo3D
  PoissonBoltzmannDemo2D
)



SET(SolverParamPath ${PACKAGE_SOURCE_DIR}/etc/SolverParameters)


ADD_TEST_BATCH(MPITests 
        COMM serial mpi
        NUM_MPI_PROCS 1
        )

PACKAGE_COPY_FILES_TO_BINARY_DIR(ExampleCopyFiles
  DEST_FILES 
                   amesos.xml 
                   anasazi-ml.xml 
                   aztec-ifpack.xml 
                   aztec-ml.xml 
                   aztec-native.xml
                   aztec.xml
                   belos-ifpack.xml 
                   belos-ml.xml 
                   bicgstab.xml
                   gmres.xml
                   nox.xml            
                   nox-amesos.xml            
                   nox-aztec.xml            
                   poissonParams.xml
                   userPrecParams.xml
  SOURCE_DIR ${SolverParamPath}
  EXEDEPS ${MPITests} ${SerialTests}
  )



