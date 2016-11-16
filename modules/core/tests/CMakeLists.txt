set( _CoreTests
  basicStatismoTest
)

foreach( _test ${_CoreTests} )
  add_executable( ${_test} ${_test}.cxx ../../ITK/MHFitting/Evaluators.h ../../ITK/MHFitting/MHFittingParameters.h ../../ITK/MHFitting/MeshOperations.h)
  target_link_libraries( ${_test} statismo_core )
endforeach()

add_test( NAME basicStatismoTest
  COMMAND basicStatismoTest
  WORKING_DIRECTORY "${TEST_WORKING_DIR}" )
