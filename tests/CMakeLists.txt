file(GLOB UNIT_TESTS_SRC
    "unit_tests/*/*.cpp"
)
file(GLOB INTEGRATION_TESTS_SRC
    "integration_tests/*/*.cpp"
)
file(GLOB REGRESSION_TESTS_SRC
    "regression_tests/*/*.cpp"
)


#set( pisces_src 
##  ../src/amoeba.c
#  ../src/C60.cpp
#  ../src/ChargeDipPol.cpp
#  ../src/ClusterAnion.cpp
#  ../src/cm_dvr.cpp
#  ../src/DPP.cpp
#  ../src/DVR.cpp
#  ../src/davdriver.cpp
#  ../src/Davidson.cpp
#  ../src/DerivDDTensor.cpp
#  ../src/GetInput.cpp
#  ../src/GeoAux.cpp
#  ../src/Gradients.cpp
#  ../src/fulldiag.cpp
#  ../src/ho_dvr.cpp
#  ../src/KE_diag.cpp
#  ../src/larnoldi.cpp
#  ../src/Model_pot.cpp
#  ../src/Molecule.cpp
#  ../src/MolPolAux.cpp
#  ../src/mtx.cpp
#  ../src/NaCl.cpp
#  ../src/optimize.cpp
#  ../src/Parameters.cpp
#  ../src/polplot.cpp
#  ../src/Potential.cpp
##  ../src/pisces.cpp
#  ../src/potfit.cpp
##  ../src/Powell.cpp
#  ../src/ReadCubeFile.cpp
#  ../src/sine_dvr.cpp
#  ../src/Small2Large.cpp
#  ../src/tsin.cpp
#  ../src/vtx_FFT.cpp
#  ../src/VectorFFTW.cpp
#  ../src/Water.cpp
#  ../src/WriteCubeFile.cpp
#  )

##############################################################################################################
# Do all the searching for pisces libraries
add_executable(PISCES_TESTS doctest_main.cpp ${UNIT_TESTS_SRC} ${INTEGRATION_TESTS_SRC} ${REGRESSION_TESTS_SRC} ) # ${pisces_src} )
target_link_libraries(PISCES_TESTS piscesf pisces_src)

if ( OPENMP_FOUND )
set_property( 
  TARGET  PISCES_TESTS 
  APPEND
  PROPERTY COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
  )
set_property( 
  TARGET  PISCES_TESTS 
  APPEND
  PROPERTY LINK_FLAGS ${OpenMP_CXX_FLAGS}
  )
else ( OPENMP_FOUND )
message(STATUS "Warning: Parallelization disabled. Compiler doesn't support OpenMP")
endif ( OPENMP_FOUND )

if(FFTW_FOUND)
  set_property(DIRECTORY PROPERTY COMPILE_DEFINITIONS HAVE_FFTW)
  include_directories(${FFTW_INCLUDES})
  target_link_libraries(PISCES_TESTS ${FFTW_LIBRARY} ${FFTW_OMP_LIBRARY})
endif(FFTW_FOUND)
if(ARPACK_FOUND)
  set_property(DIRECTORY PROPERTY COMPILE_DEFINITIONS HAVE_ARPACK)
  target_link_libraries(PISCES_TESTS ${ARPACK_LIBRARY})
endif(ARPACK_FOUND)
target_link_libraries( PISCES_TESTS ${LAPACK_LIBRARIES} )
target_link_libraries(PISCES_TESTS doctest::doctest)
target_include_directories(PISCES_TESTS PUBLIC include)
doctest_discover_tests(PISCES_TESTS)
