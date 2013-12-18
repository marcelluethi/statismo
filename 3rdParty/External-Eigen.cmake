set( Eigen3_VERSION "3.1.4" )

ExternalProject_Add( Eigen3
  URL "http://bitbucket.org/eigen/eigen/get/${Eigen3_VERSION}.tar.gz"
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/Eigen3-prefix/src/Eigen3/Eigen ${CMAKE_BINARY_DIR}/3rdParty/Eigen
)

set(EIGEN_INCLUDE_DIR ${CMAKE_BINARY_DIR}/3rdParty/Eigen/ )
