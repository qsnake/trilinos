

INCLUDE(AddTestBatch)

SET(DisabledTests
  RefinementTest3D
  FunctionalTest3D 
  transientOneZoneEddy2D
  FOSLSMagnetostatics2D
  PoissonBoltzmannDemo3D
  NoxAlgebraic 
  NonlinAlgebraic 
)

SET(SerialTests
  DiscFunc3D
  PoissonOnDisk 
  LinearPeriodic1D
  NonlinearPeriodic1D
  Poisson3D
  PoissonBoltzmannDemo1D
  BesselTest2D 
  Coupled2D
  Kepler 
  Stokes2D 
  PCDStokesCouette2D
  LinearPartialDomain 
  NonlinearPartialDomain 
  AToCDensitySample
  NonlinearPDECOTest1D
  NonlinearPDECOTest2D
  FunctionalTest 
  TwoZoneDiscFunc
  PerforatedPlate
  PerforatedPlate_QuadInt
  Stokes2D_HN_Q1
  PoissonPeano2D
  PoissonPeano3D
  StokesPeano2D 
  SecondOrderFloquet 
  MixedTest
  NavStok_Chanel_stat_HN_Nitsch
  Poisson2D_HN_Nitsch
)



SET(MPITests 
             LPTests1D
             Poisson2D
             RadDiff1D
             CrankNicolsonHeat1D
             FieldIOTest
             BlockStochPoissonTest1D
             HighOrderPoisson2D
             HighOrderPoissonBernstein2D
             HighOrderProjection2D 
             NoxPoisson1D 
             PoissonBoltzmann2D
             UserDefAlgebraic 
             TaylorHoodStokes2D
             AdvectionDiffusion2D 
             LinearVariationalTest 
             OptPoisson1D
             SimpleModelEvalTest 
             VorticityStokes2D 
             VorticityNavierStokes2D
             SpectralPoisson1D 
             SpectralSqrt 
             SteadyBurgers1D
             NitschePoisson2D 
             HelmholtzEigen2D
             MixedTest_HN
             NavStok_Chanel_stat_HN
             Poisson3DTest_HN_Q1Q1
             Poisson3DTest_HN_Q2Q1
   )



SET(SolverParamPath ${PACKAGE_SOURCE_DIR}/etc/SolverParameters)


ADD_TEST_BATCH(MPITests 
        COMM serial mpi
        )

ADD_TEST_BATCH(SerialTests 
        COMM serial mpi
        NUM_MPI_PROCS 1
        )

PACKAGE_COPY_FILES_TO_BINARY_DIR(TestCopyFiles
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
                   nox-newton.xml   
                   pcd-outer-belos.xml               
                   poissonParams.xml
                   userPrecParams.xml
  SOURCE_DIR ${SolverParamPath}
  EXEDEPS ${MPITests} ${SerialTests}
  )



# Eikonal2D needs special handling because of the fortran file

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Eikonal2D
       SOURCES Eikonal2D.cpp bessi0.f
       COMM serial mpi
       NUM_MPI_PROCS 1
   )

PACKAGE_COPY_FILES_TO_BINARY_DIR(EikonalCopyFiles
  DEST_FILES nox.xml
  SOURCE_DIR ${SolverParamPath}
  EXEDEPS Eikonal2D
  )
