INCLUDE(PackageAddExecutable)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE(
  CustomGUIExample
  SOURCES main.cpp
)

PACKAGE_COPY_FILES_TO_BINARY_DIR( CustomGUIFileCopies
  DEST_FILES
  myStyleSheet.qss myIcon.png
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/example/CustomGUIExample
  EXEDEPS CustomGUIExample
  )
