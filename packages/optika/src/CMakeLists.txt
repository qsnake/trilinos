INCLUDE(PackageLibraryMacros)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

include(${QT_USE_FILE})

SET(HDRS Optika_Dependency.hpp 
Optika_DependencySheet.hpp Optika_ArrayHelperFunctions.hpp 
Optika_GUI.hpp Optika_InvalidDependencyException.hpp
Optika_SpecificParameterEntryValidators.hpp Optika_StandardDependencies.hpp
Optika_treeitem.hpp Optika_Types.hpp Optika_ValidatorFactory.hpp Optika_Condition.hpp Optika_StandardConditions.hpp)


SET(MOC_HDRS Optika_metawindow.hpp Optika_treemodel.hpp
Optika_delegate.hpp Optika_ArrayWidget.hpp Optika_FileNameWidget.hpp
Optika_treeview.hpp)

SET(SRCS Optika_delegate.cpp Optika_Dependency.cpp 
Optika_DependencySheet.cpp Optika_ArrayHelperFunctions.cpp 
Optika_FileNameWidget.cpp Optika_GUI.cpp Optika_metawindow.cpp 
Optika_StandardDependencies.cpp Optika_treeitem.cpp 
Optika_treemodel.cpp Optika_treeview.cpp Optika_Types.cpp
Optika_ValidatorFactory.cpp Optika_SpecificParameterEntryValidators.cpp Optika_StandardConditions.cpp)
 
QT4_WRAP_CPP(MOC_SRCS ${MOC_HDRS} )


PACKAGE_ADD_LIBRARY(
	optika 
	HEADERS ${MOC_HDRS} ${HDRS} 
	SOURCES ${SRCS} ${MOC_SRCS} 
	DEPLIBS
	teuchos ${TPL_Qt_LIBRARIES}
	)

