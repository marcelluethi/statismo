if(MSVC11) #i.e. Visual Studio 2012
  # Fix for VS2012 that has _VARIADIC_MAX set to 5. Don't set too high because it increases compiler memory usage / compile-time.
  add_definitions( -D_VARIADIC_MAX=10 )
  # Fix for another VS2012 problem: not all TR1 options are automatically detected, therefore we force them here.
  add_definitions( -D BOOST_HAS_TR1 )
  add_definitions( -D BOOST_NO_0X_HDR_INITIALIZER_LIST )
endif()


set(thrift_INCLUDE_DIR "THRIFT_INCLUDE_DIR" CACHE PATH "")
set(statismoui_DIR "STATISMO_UI_DIR" CACHE PATH "")
set(thrift_LIBRARY_DIR "THRIFT_LIBRARY_DIR" CACHE PATH "")
include_directories(${statismo_INCLUDE_DIRS} ${thrift_INCLUDE_DIR} ${statismoui_DIR}/include )
link_directories( ${ITK_LIBRARY_DIRS} optimized ${statismo_LIBRARY_DIRS} ${statismoui_DIR}/lib ${thrift_LIBRARY_DIR})

set( _examples
   itkMHFittingTest
)

foreach( e ${_examples} )
  add_executable( ${e} ${e}.cxx
          Proposals.h MHFittingUtils.h)
  target_link_libraries( ${e} statismo_core ${ITK_LIBRARIES} statismo-ui thrift )
  target_link_libraries( ${e} statismo_core ${ITK_LIBRARIES} )
endforeach()
