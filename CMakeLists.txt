cmake_minimum_required(VERSION 3.16)

project(AccuSphGeom LANGUAGES C CXX)

option(ACCUSPHGEOM_PREDICATES_USE_FLOAT
       "Use float version of Shewchuk predicates"
       OFF)
option(ACCUSPHGEOM_BUILD_TESTS "Build AccuSphGeom tests" OFF)
option(ACCUSPHGEOM_ENABLE_EIGEN_ADAPTERS
       "Build Eigen adapter tests and SIMD-pack benchmarks"
       OFF)

find_path(ACCUSPHGEOM_EIGEN3_INCLUDE_DIR Eigen/Dense
  PATHS
    /usr/include/eigen3
    /usr/local/include/eigen3
)

add_library(simd_fma INTERFACE)

target_compile_options(simd_fma INTERFACE -march=native -ffp-contract=off)

if(ACCUSPHGEOM_EIGEN3_INCLUDE_DIR)
  target_include_directories(simd_fma INTERFACE ${ACCUSPHGEOM_EIGEN3_INCLUDE_DIR})
endif()

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

target_link_libraries(AccuSphGeom PUBLIC simd_fma)

if(ACCUSPHGEOM_PREDICATES_USE_FLOAT)
  target_compile_definitions(AccuSphGeom PUBLIC ACCUSPHGEOM_PREDICATES_USE_FLOAT)
endif()

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
  add_test(NAME test_gca_gca_intersection_baseline
           COMMAND test_gca_gca_intersection_baseline)

  add_executable(test_gca_constlat_intersection_baseline
    tests/test_gca_constlat_intersection_baseline.cpp
  )
  target_link_libraries(test_gca_constlat_intersection_baseline PRIVATE AccuSphGeom)
  target_compile_definitions(test_gca_constlat_intersection_baseline
    PRIVATE ACCUSPHGEOM_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
  )
  add_test(NAME test_gca_constlat_intersection_baseline
           COMMAND test_gca_constlat_intersection_baseline)

endif()

if(ACCUSPHGEOM_ENABLE_EIGEN_ADAPTERS)
  if(ACCUSPHGEOM_EIGEN3_INCLUDE_DIR)
    message(STATUS "Eigen3 found at: ${ACCUSPHGEOM_EIGEN3_INCLUDE_DIR}")

    if(ACCUSPHGEOM_BUILD_TESTS)
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
      add_test(NAME test_gca_constlat_eigen_pack
               COMMAND test_gca_constlat_eigen_pack)
    endif()

    find_package(OpenMP REQUIRED)

    add_executable(benchmark_gca_constlat_SIMDPack
      tests/performance_test/gca_constLat/benchmark_gca_constlat_SIMDPack.cpp
      tests/performance_test/gca_constLat/fp64_GCAconstLat.hh
    )

    target_link_libraries(benchmark_gca_constlat_SIMDPack
      PRIVATE
        AccuSphGeom
        OpenMP::OpenMP_CXX
    )

    target_include_directories(benchmark_gca_constlat_SIMDPack
      PRIVATE
        ${ACCUSPHGEOM_EIGEN3_INCLUDE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/performance_test/gca_constLat
    )

    if(MSVC)
      target_compile_options(benchmark_gca_constlat_SIMDPack PRIVATE /O2)
    else()
      target_compile_options(benchmark_gca_constlat_SIMDPack
        PRIVATE
          -O3
          -march=native
          -ffp-contract=off
          -fno-fast-math
      )
    endif()

    add_executable(benchmark_gca_constlat_try_SIMDPack
      tests/performance_test/gca_constLat/benchmark_gca_constlat_try_SIMDPack.cpp
      tests/performance_test/gca_constLat/fp64_GCAconstLat.hh
    )

    target_link_libraries(benchmark_gca_constlat_try_SIMDPack
      PRIVATE
        AccuSphGeom
        OpenMP::OpenMP_CXX
    )

    target_include_directories(benchmark_gca_constlat_try_SIMDPack
      PRIVATE
        ${ACCUSPHGEOM_EIGEN3_INCLUDE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/performance_test/gca_constLat
    )

    if(MSVC)
      target_compile_options(benchmark_gca_constlat_try_SIMDPack PRIVATE /O2)
    else()
      target_compile_options(benchmark_gca_constlat_try_SIMDPack
        PRIVATE
          -O3
          -march=native
          -ffp-contract=off
          -fno-fast-math
      )
    endif()

    add_executable(benchmark_gca_gca_SIMDPack
      tests/performance_test/gca_gca/benchmark_gca_gca_SIMDPack.cpp
      tests/performance_test/gca_gca/fp64_GCAGCA.hh
    )

    target_link_libraries(benchmark_gca_gca_SIMDPack
      PRIVATE
        AccuSphGeom
        OpenMP::OpenMP_CXX
    )

    target_include_directories(benchmark_gca_gca_SIMDPack
      PRIVATE
        ${ACCUSPHGEOM_EIGEN3_INCLUDE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/performance_test/gca_gca
    )

    if(MSVC)
      target_compile_options(benchmark_gca_gca_SIMDPack PRIVATE /O2)
    else()
      target_compile_options(benchmark_gca_gca_SIMDPack
        PRIVATE
          -O3
          -march=native
          -ffp-contract=off
          -fno-fast-math
      )
    endif()

    add_executable(benchmark_gca_gca_try_SIMDPack
      tests/performance_test/gca_gca/benchmark_gca_gca_try_SIMDPack.cpp
      tests/performance_test/gca_gca/fp64_GCAGCA.hh
    )

    target_link_libraries(benchmark_gca_gca_try_SIMDPack
      PRIVATE
        AccuSphGeom
        OpenMP::OpenMP_CXX
    )

    target_include_directories(benchmark_gca_gca_try_SIMDPack
      PRIVATE
        ${ACCUSPHGEOM_EIGEN3_INCLUDE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/performance_test/gca_gca
    )

    if(MSVC)
      target_compile_options(benchmark_gca_gca_try_SIMDPack PRIVATE /O2)
    else()
      target_compile_options(benchmark_gca_gca_try_SIMDPack
        PRIVATE
          -O3
          -march=native
          -ffp-contract=off
          -fno-fast-math
      )
    endif()

  else()
    message(WARNING
      "ACCUSPHGEOM_ENABLE_EIGEN_ADAPTERS=ON, but Eigen3 was not found. "
      "Eigen adapter tests and SIMD-pack benchmarks will not be built."
    )
  endif()
endif()
