# CMake tests specification 

INCLUDE(PackageAddExecutableAndTest)


IF (TPL_ENABLE_ExodusII)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
        partitionExo
        SOURCES partitionExo.cpp
        COMM serial mpi
)




PACKAGE_COPY_FILES_TO_BINARY_DIR(SundanceMeshCopyFiles
  DEST_FILES wheel.exo
  EXEDEPS partitionExo
  )

ENDIF()




