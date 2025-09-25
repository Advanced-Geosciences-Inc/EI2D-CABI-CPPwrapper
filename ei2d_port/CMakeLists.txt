cmake_minimum_required(VERSION 3.20)
project(ei2d_port LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(BUILD_SHARED_LIBS ON)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran_mods)
set(BLA_VENDOR OpenBLAS)

find_package(OpenMP)
find_package(BLAS REQUIRED)
find_package(LAPACK QUIET)

file(GLOB_RECURSE EI2D_F_SRC
  "${CMAKE_SOURCE_DIR}/apps/ei2dcore/src/*.f90"
  "${CMAKE_SOURCE_DIR}/apps/ei2dcore/src/*.f"
)

add_library(ei2d_core SHARED
  ${EI2D_F_SRC}
  apps/ei2dcore/cwrap/ei2d_c_api.f90
)

# Legacy flags + (optional) bounds checking while we stabilize
target_compile_options(ei2d_core PRIVATE
  $<$<COMPILE_LANGUAGE:Fortran>:-O3 -fopenmp -std=legacy -fallow-argument-mismatch -fcheck=bounds -g>
)
target_link_libraries(ei2d_core PUBLIC ${BLAS_LIBRARIES})
if(LAPACK_FOUND)
  target_link_libraries(ei2d_core PUBLIC ${LAPACK_LIBRARIES})
endif()
if(OpenMP_FOUND)
  target_link_libraries(ei2d_core PUBLIC OpenMP::OpenMP_Fortran)
endif()

# Allow duplicate 'CheckPoint' symbol across units
if (MINGW)
  target_link_options(ei2d_core PRIVATE -Wl,--allow-multiple-definition)
endif()

add_executable(ei2dcli apps/ei2dcli/src/ei2dcli.cpp)
target_include_directories(ei2dcli PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(ei2dcli PRIVATE ei2d_core)
if(OpenMP_FOUND)
  target_link_libraries(ei2dcli PRIVATE OpenMP::OpenMP_CXX)
endif()
