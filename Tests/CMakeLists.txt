project(statismoTests)
cmake_minimum_required(VERSION 2.6)

include_directories(
  ${STATISMO_ROOT_DIR}
  ${STATISMO_ROOT_DIR}/Representers
)

add_executable(basicStatismoTest basicStatismoTest.cpp) 
target_link_libraries(basicStatismoTest ${HDF5_LIBRARIES})

add_test(basicStatismoTest ${CMAKE_BINARY_DIR}/bin/basicStatismoTest)

add_subdirectory(modelBuilderTests)

