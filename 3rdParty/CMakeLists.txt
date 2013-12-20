project( Statismo_SuperBuild )

cmake_minimum_required( VERSION 2.8.4 )

include( ExternalProject )

set( Patches_DIR ${CMAKE_CURRENT_SOURCE_DIR}/patches )
set( Statismo_DEPENDENCIES )

# -----------------------------------------------------
# Eigen3
option( External_Eigen3 "Build Eigen3" ON )

if( ${External_Eigen3} MATCHES "ON" )
  include( External-Eigen.cmake )
  set( Statismo_DEPENDENCIES ${Statismo_DEPENDENCIES} Eigen3 )
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/3rdParty/Eigen/ DESTINATION include/Eigen )
else()
  set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../cmake )
  find_package( Eigen3 REQUIRED )
endif()

# -----------------------------------------------------
# Boost
option( External_Boost "Build Boost" ON )

if( ${External_Boost} MATCHES "ON" )
  include( External-Boost.cmake )
  set( Statismo_DEPENDENCIES ${Statismo_DEPENDENCIES} boost )
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/3rdParty/boost/ DESTINATION include/boost )
else()
  find_package( Boost REQUIRED )
endif()

# -----------------------------------------------------
# HDF5
option( External_HDF5 "Build HDF5" ON )

if( ${External_HDF5} MATCHES "ON" )
  include( External-HDF5.cmake )
  set( Statismo_DEPENDENCIES ${Statismo_DEPENDENCIES} HDF5 )

  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/3rdParty/HDF5/include/ DESTINATION include/hdf5 )
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/3rdParty/HDF5/lib/ DESTINATION lib/hdf5 )
else()
  find_package( HDF5 COMPONENTS C CXX REQUIRED )
endif()

# -----------------------------------------------------
# VTK
option( External_VTK "Build ITK" ON )

if( ${External_VTK} MATCHES "ON" )
  include( External-VTK.cmake )
  set( Statismo_DEPENDENCIES ${Statismo_DEPENDENCIES} VTK )

  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/3rdParty/VTK/include/ DESTINATION include )
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/3rdParty/VTK/bin/ DESTINATION bin )
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/3rdParty/VTK/lib/ DESTINATION lib )
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/3rdParty/VTK/share/ DESTINATION share )
else()
  find_package( VTK 6.0 REQUIRED )
  include( ${VTK_USE_FILE} )
endif()

# -----------------------------------------------------
# ITK
option( External_ITK "Build ITK" ON )

if( ${External_ITK} MATCHES "ON" )
  include( External-ITK.cmake )
  set( Statismo_DEPENDENCIES ${Statismo_DEPENDENCIES} ITK )

  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/3rdParty/ITK/include/ DESTINATION include )
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/3rdParty/ITK/bin/ DESTINATION bin )
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/3rdParty/ITK/lib/ DESTINATION lib )
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/3rdParty/ITK/share/ DESTINATION share )
else()
  find_package( ITK 4.5 REQUIRED )
  include( ${ITK_USE_FILE} )
endif()

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
    -DBoost_INCLUDE_DIR:PATH=${Boost_INCLUDE_DIR}
    -DITK_DIR:PATH=${ITK_DIR}
    -DVTK_DIR:PATH=${VTK_DIR}
    -DHDF5_DIR:PATH=${HDF5_DIR}
    -DEIGEN3_INCLUDE_DIR:PATH=${EIGEN3_INCLUDE_DIR}
    -DBUILD_ITK_EXAMPLES:BOOL=ON
    -DBUILD_VTK_EXAMPLES:BOOL=ON
    -DBUILD_STATISMO_TESTS:BOOL=ON
    -DBUILD_REPRESENTER_TESTS:BOOL=ON
  INSTALL_COMMAND ""
)
