cmake_minimum_required(VERSION 3.16)

project(AccuSphGeom LANGUAGES C CXX)

option(ACCUSPHGEOM_BUILD_TESTS "Build AccuSphGeom tests" ON)

add_library(AccuSphGeom STATIC
  third_party/predicates.c
  third_party/geogram_psm/MultiPrecision_psm.cpp
  src/algorithms/point_in_polygon_sphere.cpp
  src/predicates/quadruple3d_exact.cpp
)

target_compile_features(AccuSphGeom PUBLIC cxx_std_17)

target_include_directories(AccuSphGeom
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/geogram_psm>
    $<INSTALL_INTERFACE:include>
)

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

  foreach(test_name IN ITEMS
      test_library_usage
      test_quadruple3d_basic
      test_pip_robust
      test_pip_complicated
      test_constructions)
    add_executable(${test_name} tests/${test_name}.cpp)
    target_link_libraries(${test_name} PRIVATE AccuSphGeom)
    add_test(NAME ${test_name} COMMAND ${test_name})
  endforeach()
endif()
