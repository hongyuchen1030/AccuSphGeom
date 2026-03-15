#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "spip/core/types.hpp"
#include "spip/kernels/pip_kernel_eft.hpp"

namespace spip::pip {

enum class Location : std::uint8_t { Outside, Inside, OnVertex, OnEdge };

// No global IDs: use the existing internal degeneracy rule.
Location point_in_polygon_sphere(const double* q,
                                 const double* const* poly,
                                 std::size_t n);

// With global IDs: use the global-ID-based SoS perturbation rule.
Location point_in_polygon_sphere(const double* q,
                                 const double* const* poly,
                                 const std::int64_t* global_vertex_ids,
                                 std::size_t n);

Location point_in_polygon_sphere(const std::array<double, 3>& q,
                                 const std::vector<std::array<double, 3>>& poly);

Location point_in_polygon_sphere(
    const std::array<double, 3>& q,
    const std::vector<std::array<double, 3>>& poly,
    const std::vector<std::int64_t>& global_vertex_ids);

Location point_in_polygon_sphere(const std::vector<double>& q,
                                 const std::vector<std::vector<double>>& poly);

Location point_in_polygon_sphere(
    const std::vector<double>& q,
    const std::vector<std::vector<double>>& poly,
    const std::vector<std::int64_t>& global_vertex_ids);

// EFT overloads. These are templated and therefore implemented in this header.
template <typename T>
Location point_in_polygon_sphere(const V3_T<T>& q,
                                 const V3_T<T>* poly,
                                 std::size_t n);

template <typename T>
inline Location point_in_polygon_sphere(const V3_T<T>& q,
                                        const std::vector<V3_T<T>>& poly) {
  return point_in_polygon_sphere(q, poly.data(), poly.size());
}

namespace detail {

template <typename T>
inline void require_valid_polygon(const V3_T<T>* poly, std::size_t n) {
  if (poly == nullptr) {
    throw std::invalid_argument("point_in_polygon_sphere: poly is null");
  }
  if (n < 3) {
    throw std::invalid_argument(
        "point_in_polygon_sphere: polygon must have >= 3 vertices");
  }
}

template <typename T>
inline bool equal3(const V3_T<T>& a, const V3_T<T>& b) {
  return (a[0] == b[0]) && (a[1] == b[1]) && (a[2] == b[2]);
}

// Return true iff q lies on the closed non-antipodal minor arc AB.
//
// The test consists of two conditions:
//
// (1) q lies on the supporting great circle of AB:
//       orient(A, B, q) = 0
//
// (2) q lies between A and B on the minor arc. Equivalently, the normals
//     A x q and q x B do not point in opposite directions:
//       (A x q) . (q x B) >= 0
//
// The second condition is evaluated as the scalar quadruple product
// quadruple3d(A, q, q, B).
template <typename T>
inline bool on_minor_arc_eft(const V3_T<T>& q,
                             const V3_T<T>& A,
                             const V3_T<T>& B) {
  using Kernel = spip::kernels::PIPKernelEFT;
  using Sign = spip::predicates::Sign;

  if (equal3(A, B)) {
    return false;
  }
  if (Kernel::orient3d_on_sphere(A, B, q) != Sign::Zero) {
    return false;
  }
  return Kernel::quadruple3d(A, q, q, B) != Sign::Negative;
}

template <typename T>
inline bool point_on_polygon_edge_eft(const V3_T<T>& q,
                                      const V3_T<T>* poly,
                                      std::size_t n) {
  const V3_T<T>* A = &poly[n - 1];
  for (std::size_t i = 0; i < n; ++i) {
    const V3_T<T>* B = &poly[i];
    if (on_minor_arc_eft(q, *A, *B)) {
      return true;
    }
    A = B;
  }
  return false;
}

// Construct a deterministic perturbation of the antipode of q.
//
// This endpoint is used only by the non-SoS EFT overload. A small perturbation
// is added to the antipode and the result is renormalized to the unit sphere.
template <typename T>
inline V3_T<T> make_perturbed_antipode_eft(const V3_T<T>& q) {
  V3_T<T> r = -q;

  const T ax = std::abs(r[0]);
  const T ay = std::abs(r[1]);
  const T az = std::abs(r[2]);

  const T eps = T(1e-15);
  if (ax <= ay && ax <= az) {
    r[0] += eps;
  } else if (ay <= ax && ay <= az) {
    r[1] += eps;
  } else {
    r[2] += eps;
  }

  const T n2 = r.dot(r);
  if (n2 > T(0)) {
    r /= std::sqrt(n2);
  }
  return r;
}

// Return true iff the polygon edge AB contributes one crossing to the ray
// crossing count for the minor arc qR, where R is a perturbed antipode.
//
// The decision follows the half-open convention of Hormann and Agathos (2001).
template <typename T>
inline bool counts_as_ray_crossing_half_open_eft(const V3_T<T>& A,
                                                 const V3_T<T>& B,
                                                 const V3_T<T>& q,
                                                 const V3_T<T>& R) {
  using Kernel = spip::kernels::PIPKernelEFT;
  using Sign = spip::predicates::Sign;

  const Sign s_qR_A = Kernel::orient3d_on_sphere(q, R, A);
  const Sign s_qR_B = Kernel::orient3d_on_sphere(q, R, B);
  const Sign s_AB_q = Kernel::orient3d_on_sphere(A, B, q);
  const Sign s_AB_R = Kernel::orient3d_on_sphere(A, B, R);

  const bool below_a = (s_qR_A == Sign::Negative);
  const bool below_b = (s_qR_B == Sign::Negative);
  if (!(below_a ^ below_b)) {
    return false;
  }

  if (s_AB_q == Sign::Zero) {
    throw std::domain_error(
        "counts_as_ray_crossing_half_open_eft: q lies on AB great circle "
        "(boundary case)");
  }
  if (s_AB_R == Sign::Zero) {
    throw std::domain_error(
        "counts_as_ray_crossing_half_open_eft: R lies on AB great circle "
        "(ray degeneracy)");
  }

  const bool q_side = (s_AB_q == Sign::Negative);
  const bool r_side = (s_AB_R == Sign::Negative);
  return (q_side ^ r_side);
}

}  // namespace detail

template <typename T>
inline Location point_in_polygon_sphere(const V3_T<T>& q,
                                        const V3_T<T>* poly,
                                        std::size_t n) {
  detail::require_valid_polygon(poly, n);

  for (std::size_t i = 0; i < n; ++i) {
    if (detail::equal3(q, poly[i])) {
      return Location::OnVertex;
    }
  }

  if (detail::point_on_polygon_edge_eft(q, poly, n)) {
    return Location::OnEdge;
  }

  const V3_T<T> R = detail::make_perturbed_antipode_eft(q);

  bool inside = false;
  const V3_T<T>* A = &poly[n - 1];
  for (std::size_t i = 0; i < n; ++i) {
    const V3_T<T>* B = &poly[i];
    if (detail::counts_as_ray_crossing_half_open_eft(*A, *B, q, R)) {
      inside = !inside;
    }
    A = B;
  }

  return inside ? Location::Inside : Location::Outside;
}

}  // namespace spip::pip