#pragma once

#include <array>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>

#include "accusphgeom/algorithms/point_in_polygon_sphere.hpp"
#include "accusphgeom/numeric/constants.hpp"
#include "accusphgeom/numeric/eft.hpp"

namespace accusphgeom::algorithms {

enum class FaceLocationLabel : std::uint8_t {
  Local,
  NorthPolarCapCandidate,
  SouthPolarCapCandidate,
};

inline constexpr double default_polar_cap_lat_deg = 80.0;
inline constexpr double default_endpoint_lat_snap_tol_deg = 0.0001;

template <typename T>
using MaskType = decltype(std::declval<T>() < std::declval<T>());

template <typename T>
struct FaceLocationInfo {
  T face_z_max{};
  T face_z_min{};
  MaskType<T> north_pole_candidate_mask{};
  MaskType<T> south_pole_candidate_mask{};
  MaskType<T> local_mask{};
};

template <typename T>
struct LatLonBounds {
  T lat_min{};
  T lat_max{};
  T lon_min{};
  T lon_max{};
  bool longitude_wraps = false;
};

template <typename T, std::size_t N>
inline FaceLocationInfo<T> get_face_location_info(
    const std::array<numeric::Vec3<T>, N>& face_vertices,
    T polar_cap_lat_deg = static_cast<T>(default_polar_cap_lat_deg)) {
  using std::max;
  using std::min;
  using std::sin;
  using std::sqrt;

  T face_z_max = -std::numeric_limits<T>::infinity();
  T face_z_min =  std::numeric_limits<T>::infinity();

  for (std::size_t i = 0; i < N; ++i) {
    const numeric::Vec3<T>& x1 = face_vertices[i];
    const numeric::Vec3<T>& x2 = face_vertices[(i + 1) % N];

    const T z1 = x1[2];
    const T z2 = x2[2];
    const T d = x1[0] * x2[0] + x1[1] * x2[1] + x1[2] * x2[2];

    const T a_raw = (z1 * d - z2) / ((z1 + z2) * (d - T(1)));
    const T a = min(max(a_raw, T(0)), T(1));

    const T one_minus_a = T(1) - a;
    const T y0 = one_minus_a * x1[0] + a * x2[0];
    const T y1 = one_minus_a * x1[1] + a * x2[1];
    const T y2 = one_minus_a * x1[2] + a * x2[2];
    const T norm = sqrt(y0 * y0 + y1 * y1 + y2 * y2);
    const T z_ext = y2 / norm;

    const T z_edge_max = max(z1, z2);
    const T z_edge_min = min(z1, z2);

    const MaskType<T> use_ext = (a_raw > T(0)) & (a_raw < T(1));

    const T z_max = use_ext ? z_ext : z_edge_max;
    const T z_min = use_ext ? z_ext : z_edge_min;

    face_z_max = max(face_z_max, z_max);
    face_z_min = min(face_z_min, z_min);
  }

  const T polar_cap_z =
      sin(polar_cap_lat_deg * numeric::pi<T> / T(180));

  const MaskType<T> north_pole_candidate_mask = (face_z_max >= polar_cap_z);
  const MaskType<T> south_pole_candidate_mask = (face_z_min <= -polar_cap_z);
  const MaskType<T> local_mask = !(north_pole_candidate_mask | south_pole_candidate_mask);

  return {
      face_z_max,
      face_z_min,
      north_pole_candidate_mask,
      south_pole_candidate_mask,
      local_mask
  };
}

template <typename T, std::size_t N>
inline LatLonBounds<T> generate_lat_lon_bounds_local(
    const std::array<numeric::Vec3<T>, N>& face_vertices,
    const FaceLocationInfo<T>& info,
    T endpoint_lat_snap_tol_deg =
        static_cast<T>(default_endpoint_lat_snap_tol_deg)) {
  using std::abs;
  using std::asin;
  using std::atan2;
  using std::max;
  using std::min;

  LatLonBounds<T> bounds{};

  const T rad_to_deg = T(180) / numeric::pi<T>;

  T endpoint_lat_max = -std::numeric_limits<T>::infinity();
  T endpoint_lat_min =  std::numeric_limits<T>::infinity();
  T endpoint_lon_max = -std::numeric_limits<T>::infinity();
  T endpoint_lon_min =  std::numeric_limits<T>::infinity();

  for (std::size_t i = 0; i < N; ++i) {
    const numeric::Vec3<T>& x = face_vertices[i];

    const T lat_deg = asin(x[2]) * rad_to_deg;
    const T lon_deg = atan2(x[1], x[0]) * rad_to_deg;

    endpoint_lat_max = max(endpoint_lat_max, lat_deg);
    endpoint_lat_min = min(endpoint_lat_min, lat_deg);
    endpoint_lon_max = max(endpoint_lon_max, lon_deg);
    endpoint_lon_min = min(endpoint_lon_min, lon_deg);
  }

  T lat_max = asin(info.face_z_max) * rad_to_deg;
  T lat_min = asin(info.face_z_min) * rad_to_deg;

  const MaskType<T> snap_lat_max_mask =
      (abs(lat_max - endpoint_lat_max) <= endpoint_lat_snap_tol_deg);
  const MaskType<T> snap_lat_min_mask =
      (abs(lat_min - endpoint_lat_min) <= endpoint_lat_snap_tol_deg);

  lat_max = snap_lat_max_mask ? endpoint_lat_max : lat_max;
  lat_min = snap_lat_min_mask ? endpoint_lat_min : lat_min;

  bounds.lat_max = lat_max;
  bounds.lat_min = lat_min;
  bounds.lon_max = endpoint_lon_max;
  bounds.lon_min = endpoint_lon_min;

  return bounds;
}

template <typename T, std::size_t N, typename GlobalId = std::int64_t>
inline LatLonBounds<T> generate_lat_lon_bounds_pole(
    const std::array<numeric::Vec3<T>, N>& face_vertices,
    const FaceLocationInfo<T>& info,
    const GlobalId* global_vertex_ids = nullptr,
    const GlobalId* pole_query_id = nullptr,
    const GlobalId* waypoint_id = nullptr,
    T endpoint_lat_snap_tol_deg =
        static_cast<T>(default_endpoint_lat_snap_tol_deg)) {
  static_assert(std::is_integral_v<GlobalId>,
                "GlobalId must be an integral type");

  LatLonBounds<T> bounds{};

  const T rad_to_deg = T(180) / numeric::pi<T>;

  constexpr GlobalId kNorthPoleQueryId =
      static_cast<GlobalId>(std::numeric_limits<GlobalId>::max() -
                            static_cast<GlobalId>(4));
  constexpr GlobalId kSouthPoleQueryId =
      static_cast<GlobalId>(std::numeric_limits<GlobalId>::max() -
                            static_cast<GlobalId>(3));
  constexpr GlobalId kWaypointId =
      static_cast<GlobalId>(std::numeric_limits<GlobalId>::max() -
                            static_cast<GlobalId>(2));

  constexpr double north_pole[3] = {0.0, 0.0, 1.0};
  constexpr double south_pole[3] = {0.0, 0.0, -1.0};

  std::array<const double*, N> face_ptrs{};
  for (std::size_t i = 0; i < N; ++i) {
    face_ptrs[i] = face_vertices[i].data();
  }

  std::array<GlobalId, N> local_ids{};
  const GlobalId* ids = global_vertex_ids;
  if (ids == nullptr) {
    for (std::size_t i = 0; i < N; ++i) {
      local_ids[i] = static_cast<GlobalId>(i);
    }
    ids = local_ids.data();
  }

  const GlobalId north_q_id =
      pole_query_id ? *pole_query_id : kNorthPoleQueryId;
  const GlobalId south_q_id =
      pole_query_id ? *pole_query_id : kSouthPoleQueryId;

  const GlobalId r_id = waypoint_id ? *waypoint_id : kWaypointId;
  (void)r_id;

  const Location north_loc = point_in_polygon_sphere(
      north_pole, north_q_id, face_ptrs.data(), ids, N);
  const Location south_loc = point_in_polygon_sphere(
      south_pole, south_q_id, face_ptrs.data(), ids, N);

  if (north_loc == Location::Outside && south_loc == Location::Outside) {
    return generate_lat_lon_bounds_local(face_vertices, info,
                                         endpoint_lat_snap_tol_deg);
  }

  T endpoint_lat_max = -std::numeric_limits<T>::infinity();
  T endpoint_lat_min = std::numeric_limits<T>::infinity();
  T endpoint_lon_max = -std::numeric_limits<T>::infinity();
  T endpoint_lon_min = std::numeric_limits<T>::infinity();

  for (std::size_t i = 0; i < N; ++i) {
    const auto& x = face_vertices[i];
    const T lat_deg =
        std::asin(std::clamp(x[2], T(-1), T(1))) * rad_to_deg;
    const T lon_deg = std::atan2(x[1], x[0]) * rad_to_deg;

    endpoint_lat_max = std::max(endpoint_lat_max, lat_deg);
    endpoint_lat_min = std::min(endpoint_lat_min, lat_deg);
    endpoint_lon_max = std::max(endpoint_lon_max, lon_deg);
    endpoint_lon_min = std::min(endpoint_lon_min, lon_deg);
  }

  T lat_max = std::asin(std::clamp(info.face_z_max, T(-1), T(1))) * rad_to_deg;
  T lat_min = std::asin(std::clamp(info.face_z_min, T(-1), T(1))) * rad_to_deg;

  if (std::abs(lat_max - endpoint_lat_max) <= endpoint_lat_snap_tol_deg) {
    lat_max = endpoint_lat_max;
  }
  if (std::abs(lat_min - endpoint_lat_min) <= endpoint_lat_snap_tol_deg) {
    lat_min = endpoint_lat_min;
  }

  if (north_loc != Location::Outside) {
    bounds.lat_min = lat_min;
    bounds.lat_max = T(90);
    bounds.lon_min = (north_loc == Location::Inside) ? T(0) : endpoint_lon_min;
    bounds.lon_max =
        (north_loc == Location::Inside) ? T(360) : endpoint_lon_max;
    bounds.longitude_wraps = (north_loc == Location::Inside);
    return bounds;
  }

  bounds.lat_min = T(-90);
  bounds.lat_max = lat_max;
  bounds.lon_min = (south_loc == Location::Inside) ? T(0) : endpoint_lon_min;
  bounds.lon_max =
      (south_loc == Location::Inside) ? T(360) : endpoint_lon_max;
  bounds.longitude_wraps = (south_loc == Location::Inside);
  return bounds;
}

}  // namespace accusphgeom::algorithms
