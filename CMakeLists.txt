cmake_minimum_required(VERSION 3.16)

project(AccuSphGeom LANGUAGES C CXX)

option(ACCUSPHGEOM_PREDICATES_USE_FLOAT
       "Use float version of Shewchuk predicates"
       OFF)
option(ACCUSPHGEOM_BUILD_TESTS "Build AccuSphGeom tests" ON)

find_path(ACCUSPHGEOM_EIGEN3_INCLUDE_DIR Eigen/Dense
  PATHS
    /usr/include/eigen3
    /usr/local/include/eigen3
)

add_library(simd_fma INTERFACE)
# Compilation options for FMA instructions
target_compile_options(simd_fma INTERFACE -march=native -ffp-contract=off)
# A note on -ffp-contract:
#   The default on GCC (at least with -std=gnu17 instead of -std=c++17) is
#   -fp-contract=fast, which can contract floating point multiplications and
#   additions spanning across **multiple statements** into FMA instructions.
#   This can be surprising and change the roundoff behavior. This is especially
#   unsafe for EFT algorithms, which rely fundamentally on the IEEE 754
#   roundoff behavior implied by the sequence of statements as literally
#   written in the source file.

add_library(AccuSphGeom STATIC
  third_party/predicates.c
  third_party/geogram_psm/MultiPrecision_psm.cpp
  src/algorithms/point_in_polygon_sphere.cpp
  src/predicates/quadruple3d_exact.cpp
)

target_compile_features(AccuSphGeom PUBLIC cxx_std_17)

# Public headers live in include/, and orient3d.hpp includes "predicates.h",
# so we expose third_party/ as an include dir for consumers too.
target_include_directories(AccuSphGeom
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/geogram_psm>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(AccuSphGeom PUBLIC simd_fma)

if(ACCUSPHGEOM_PREDICATES_USE_FLOAT)
  target_compile_definitions(AccuSphGeom PUBLIC ACCUSPHGEOM_PREDICATES_USE_FLOAT)
endif()

# Reasonable warnings and strict floating-point semantics.
if(MSVC)
  target_compile_options(AccuSphGeom PRIVATE /W4)
else()
  target_compile_options(AccuSphGeom
    PRIVATE
      -Wall
      -Wextra
      -Wpedantic
      -O3
      -ffp-contract=off
      -fno-fast-math
  )
endif()

if(ACCUSPHGEOM_BUILD_TESTS)
  enable_testing()

  add_executable(test_library_usage
    tests/test_library_usage.cpp
  )
  target_link_libraries(test_library_usage PRIVATE AccuSphGeom)
  add_test(NAME test_library_usage COMMAND test_library_usage)

  add_executable(test_psm_compile
    tests/test_psm_compile.cpp
  )
  target_link_libraries(test_psm_compile PRIVATE AccuSphGeom)
  add_test(NAME test_psm_compile COMMAND test_psm_compile)

  add_executable(test_quadruple3d_basic
    tests/test_quadruple3d_basic.cpp
  )
  target_link_libraries(test_quadruple3d_basic PRIVATE AccuSphGeom)
  add_test(NAME test_quadruple3d_basic COMMAND test_quadruple3d_basic)

  add_executable(test_pip_robust
    tests/test_pip_robust.cpp
  )
  target_link_libraries(test_pip_robust PRIVATE AccuSphGeom)
  add_test(NAME test_pip_robust COMMAND test_pip_robust)

  add_executable(test_pip_complicated
    tests/test_pip_complicated.cpp
  )
  target_link_libraries(test_pip_complicated PRIVATE AccuSphGeom)
  add_test(NAME test_pip_complicated COMMAND test_pip_complicated)

  add_executable(test_gca_gca_intersection_baseline
    tests/test_gca_gca_intersection_baseline.cpp
  )
  target_link_libraries(test_gca_gca_intersection_baseline PRIVATE AccuSphGeom)
  target_compile_definitions(test_gca_gca_intersection_baseline
    PRIVATE ACCUSPHGEOM_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
  )
  add_test(NAME test_gca_gca_intersection_baseline COMMAND test_gca_gca_intersection_baseline)

  add_executable(test_gca_constlat_intersection_baseline
    tests/test_gca_constlat_intersection_baseline.cpp
  )
  target_link_libraries(test_gca_constlat_intersection_baseline PRIVATE AccuSphGeom)
  target_compile_definitions(test_gca_constlat_intersection_baseline
    PRIVATE ACCUSPHGEOM_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
  )
  add_test(NAME test_gca_constlat_intersection_baseline COMMAND test_gca_constlat_intersection_baseline)

  if(ACCUSPHGEOM_EIGEN3_INCLUDE_DIR)
    add_executable(test_gca_gca_eigen_pack
      tests/test_gca_gca_eigen_pack.cpp
    )
    target_link_libraries(test_gca_gca_eigen_pack PRIVATE AccuSphGeom)
    target_include_directories(test_gca_gca_eigen_pack
      PRIVATE ${ACCUSPHGEOM_EIGEN3_INCLUDE_DIR}
    )
    add_test(NAME test_gca_gca_eigen_pack COMMAND test_gca_gca_eigen_pack)

    add_executable(test_gca_constlat_eigen_pack
      tests/test_gca_constlat_eigen_pack.cpp
    )
    target_link_libraries(test_gca_constlat_eigen_pack PRIVATE AccuSphGeom)
    target_include_directories(test_gca_constlat_eigen_pack
      PRIVATE ${ACCUSPHGEOM_EIGEN3_INCLUDE_DIR}
    )
    add_test(NAME test_gca_constlat_eigen_pack COMMAND test_gca_constlat_eigen_pack)
  endif()

endif()
