find_package(Doxygen REQUIRED)
set( SRCROOT ${PROJECT_SOURCE_DIR}/src )
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY ) 
add_custom_target(doc ${DOXYGEN_EXECUTABLE})
