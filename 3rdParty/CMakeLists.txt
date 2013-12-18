project( Statismo_SuperBuild )

cmake_minimum_required( VERSION 2.8.4 )

include( ExternalProject )

set( Patches_DIR ${CMAKE_CURRENT_SOURCE_DIR}/patches )

include( External-Eigen.cmake )
include( External-Boost.cmake )
include( External-HDF5.cmake )
include( External-ITK.cmake )
include( External-VTK.cmake )

set( Statismo_DEPENDENCIES Eigen3 boost HDF5 ITK VTK )

ExternalProject_Add( Statismo
  DEPENDS ${Statismo_DEPENDENCIES}
  DOWNLOAD_COMMAND ""
  UPDATE_COMMAND ""
  PATCH_COMMAND ""
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..
  BINARY_DIR Statismo-build
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    ${ep_common_args}
    -DITK_DIR:PATH=${ITK_DIR}
    -DVTK_DIR:PATH=${VTK_DIR}
    -DBUILD_ITK_EXAMPLES:BOOL=ON
    -DBUILD_VTK_EXAMPLES:BOOL=ON
    -DBUILD_STATISMO_TESTS:BOOL=ON
    -DBUILD_REPRESENTER_TESTS:BOOL=ON
  INSTALL_COMMAND ""
)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/3rdParty/HDF5/include/ DESTINATION include/hdf5 FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/3rdParty/HDF5/include/cpp/ DESTINATION include/hdf5/cpp FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/3rdParty/HDF5/lib/ DESTINATION lib/hdf5 FILES_MATCHING PATTERN "*hdf5*")
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/3rdParty/HDF5/bin/ DESTINATION bin FILES_MATCHING PATTERN "*hdf5*")

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/3rdParty/boost/ DESTINATION include/boost )

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/3rdParty/Eigen/ DESTINATION include/Eigen )
