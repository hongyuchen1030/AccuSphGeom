#include "fp64_GCAconstLat.hh"

#include "accusphgeom/adapters/eigen/numeric.hpp"

namespace accusphgeom::performance_test {

template <typename T>
static accusphgeom::numeric::Vec3<T> simd_cross(
    const accusphgeom::numeric::Vec3<T>& v1,
    const accusphgeom::numeric::Vec3<T>& v2) {
  accusphgeom::numeric::Vec3<T> n;
  n[0] = v1[1] * v2[2] - v1[2] * v2[1];  // x = Ay * Bz - Az * By
  n[1] = v1[2] * v2[0] - v1[0] * v2[2];  // y = Az * Bx - Ax * Bz
  n[2] = v1[0] * v2[1] - v1[1] * v2[0];  // z = Ax * By - Ay * Bx
  return n;
}

template <typename T>
std::tuple<T, T> fp64_gca_constlat_pure_fp_xy(
    const accusphgeom::numeric::Vec3<T>& pointA,
    const accusphgeom::numeric::Vec3<T>& pointB,
    const T& constZ) {
  // Calculate the cross product, which gives us the vector n
  accusphgeom::numeric::Vec3<T> n = simd_cross(pointA, pointB);

  // Extract nx, ny, nz components from vector n
  T nx = n[0];
  T ny = n[1];
  T nz = n[2];

  // Calculate s (as 's' in the formula)
  T nx_squared = nx * nx;
  T ny_squared = ny * ny;
  T nz_squared = nz * nz;
  T nx_squared_plus_ny_squared = nx_squared + ny_squared;
  T norm_n_squared = nx_squared_plus_ny_squared + nz_squared;  // ||n||^2
  T s_tilde =
      sqrt(nx_squared_plus_ny_squared - norm_n_squared * constZ * constZ);

  // Calculate p_x and p_y, which are the x and y coordinates of the
  // intersection.
  T p_x = -(constZ * nx * nz + s_tilde * ny) / nx_squared_plus_ny_squared;

  T p_y =
      -(constZ * ny * nz - s_tilde * nx) / nx_squared_plus_ny_squared;

  return {p_x, p_y};
}

template std::tuple<double, double> fp64_gca_constlat_pure_fp_xy<double>(
    const accusphgeom::numeric::Vec3<double>& pointA,
    const accusphgeom::numeric::Vec3<double>& pointB,
    const double& constZ);

template std::tuple<EigenPack<2>, EigenPack<2>>
fp64_gca_constlat_pure_fp_xy<EigenPack<2>>(
    const accusphgeom::numeric::Vec3<EigenPack<2>>& pointA,
    const accusphgeom::numeric::Vec3<EigenPack<2>>& pointB,
    const EigenPack<2>& constZ);

template std::tuple<EigenPack<4>, EigenPack<4>>
fp64_gca_constlat_pure_fp_xy<EigenPack<4>>(
    const accusphgeom::numeric::Vec3<EigenPack<4>>& pointA,
    const accusphgeom::numeric::Vec3<EigenPack<4>>& pointB,
    const EigenPack<4>& constZ);

template std::tuple<EigenPack<8>, EigenPack<8>>
fp64_gca_constlat_pure_fp_xy<EigenPack<8>>(
    const accusphgeom::numeric::Vec3<EigenPack<8>>& pointA,
    const accusphgeom::numeric::Vec3<EigenPack<8>>& pointB,
    const EigenPack<8>& constZ);

}  // namespace accusphgeom::performance_test
