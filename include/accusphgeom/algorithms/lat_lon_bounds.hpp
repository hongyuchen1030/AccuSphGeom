#pragma once

#include <array>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

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
struct FaceLocationInfo {
  FaceLocationLabel label = FaceLocationLabel::Local;
  T face_z_max{};
  T face_z_min{};
};

template <typename T>
struct LatLonBounds {
  T lat_min{};
  T lat_max{};
  T lon_min{};
  T lon_max{};
  bool longitude_wraps = false;
};

inline FaceLocationInfo<double> get_face_location_info(
    const double* const* face_vertices,
    std::size_t n,
    double polar_cap_lat_deg = default_polar_cap_lat_deg) {
  double face_z_max = -std::numeric_limits<double>::infinity();
  double face_z_min = std::numeric_limits<double>::infinity();

  for (std::size_t i = 0; i < n; ++i) {
    const double* x1 = face_vertices[i];
    const double* x2 = face_vertices[(i + 1) % n];

    const double z1 = x1[2];
    const double z2 = x2[2];
    const double d = x1[0] * x2[0] + x1[1] * x2[1] + x1[2] * x2[2];
    double a = (z1 * d - z2) / ((z1 + z2) * (d - 1.0));
    a = std::clamp(a, 0.0, 1.0);

    const double one_minus_a = 1.0 - a;
    const double y0 = one_minus_a * x1[0] + a * x2[0];
    const double y1 = one_minus_a * x1[1] + a * x2[1];
    const double y2 = one_minus_a * x1[2] + a * x2[2];
    const double norm = std::sqrt(y0 * y0 + y1 * y1 + y2 * y2);
    const double za = y2 / norm;

    face_z_max = std::max(face_z_max, za);
    face_z_min = std::min(face_z_min, za);
  }

  const double polar_cap_z =
      std::sin(polar_cap_lat_deg * numeric::pi<double> / 180.0);

  FaceLocationLabel label = FaceLocationLabel::Local;
  if (face_z_max >= polar_cap_z) {
    label = FaceLocationLabel::NorthPolarCapCandidate;
  } else if (face_z_min <= -polar_cap_z) {
    label = FaceLocationLabel::SouthPolarCapCandidate;
  }

  return {label, face_z_max, face_z_min};
}

inline FaceLocationLabel get_face_location_label(
    const double* const* face_vertices,
    std::size_t n,
    double polar_cap_lat_deg = default_polar_cap_lat_deg) {
  return get_face_location_info(face_vertices, n, polar_cap_lat_deg).label;
}

inline LatLonBounds<double> generate_lat_lon_bounds_local(
    const double* const* face_vertices,
    std::size_t n,
    double endpoint_lat_snap_tol_deg = default_endpoint_lat_snap_tol_deg) {
  const FaceLocationInfo<double> info = get_face_location_info(face_vertices, n);
  LatLonBounds<double> bounds{};

  const double rad_to_deg = 180.0 / numeric::pi<double>;

  double endpoint_lat_max = -std::numeric_limits<double>::infinity();
  double endpoint_lat_min = std::numeric_limits<double>::infinity();
  double endpoint_lon_max = -std::numeric_limits<double>::infinity();
  double endpoint_lon_min = std::numeric_limits<double>::infinity();

  for (std::size_t i = 0; i < n; ++i) {
    const double* x = face_vertices[i];
    const double z = std::clamp(x[2], -1.0, 1.0);
    const double lat_deg = std::asin(z) * rad_to_deg;
    const double lon_deg = std::atan2(x[1], x[0]) * rad_to_deg;

    endpoint_lat_max = std::max(endpoint_lat_max, lat_deg);
    endpoint_lat_min = std::min(endpoint_lat_min, lat_deg);
    endpoint_lon_max = std::max(endpoint_lon_max, lon_deg);
    endpoint_lon_min = std::min(endpoint_lon_min, lon_deg);
  }

  double lat_max =
      std::asin(std::clamp(info.face_z_max, -1.0, 1.0)) * rad_to_deg;
  double lat_min =
      std::asin(std::clamp(info.face_z_min, -1.0, 1.0)) * rad_to_deg;

  if (std::abs(lat_max - endpoint_lat_max) <= endpoint_lat_snap_tol_deg) {
    lat_max = endpoint_lat_max;
  }
  if (std::abs(lat_min - endpoint_lat_min) <= endpoint_lat_snap_tol_deg) {
    lat_min = endpoint_lat_min;
  }

  bounds.lat_max = lat_max;
  bounds.lat_min = lat_min;
  bounds.lon_max = endpoint_lon_max;
  bounds.lon_min = endpoint_lon_min;
  bounds.longitude_wraps = false;
  return bounds;
}

inline LatLonBounds<double> generate_lat_lon_bounds_pole(
    const double* const* face_vertices,
    std::size_t n,
    double polar_cap_lat_deg = default_polar_cap_lat_deg,
    double endpoint_lat_snap_tol_deg = default_endpoint_lat_snap_tol_deg) {
  const FaceLocationInfo<double> info =
      get_face_location_info(face_vertices, n, polar_cap_lat_deg);
  LatLonBounds<double> bounds{};

  constexpr std::int64_t kNorthPoleQueryId =
      std::numeric_limits<std::int64_t>::max() - 4;
  constexpr std::int64_t kSouthPoleQueryId =
      std::numeric_limits<std::int64_t>::max() - 3;

  constexpr double north_pole[3] = {0.0, 0.0, 1.0};
  constexpr double south_pole[3] = {0.0, 0.0, -1.0};

  std::vector<std::int64_t> local_ids;
  local_ids.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    local_ids.push_back(static_cast<std::int64_t>(i));
  }

  const Location north_loc = point_in_polygon_sphere(
      north_pole, kNorthPoleQueryId, face_vertices, local_ids.data(), n);
  const Location south_loc = point_in_polygon_sphere(
      south_pole, kSouthPoleQueryId, face_vertices, local_ids.data(), n);

  if (north_loc == Location::Outside && south_loc == Location::Outside) {
    return generate_lat_lon_bounds_local(face_vertices, n,
                                         endpoint_lat_snap_tol_deg);
  }

  const double rad_to_deg = 180.0 / numeric::pi<double>;

  double endpoint_lat_max = -std::numeric_limits<double>::infinity();
  double endpoint_lat_min = std::numeric_limits<double>::infinity();
  double endpoint_lon_max = -std::numeric_limits<double>::infinity();
  double endpoint_lon_min = std::numeric_limits<double>::infinity();

  for (std::size_t i = 0; i < n; ++i) {
    const double* x = face_vertices[i];
    const double lat_deg =
        std::asin(std::clamp(x[2], -1.0, 1.0)) * rad_to_deg;
    const double lon_deg = std::atan2(x[1], x[0]) * rad_to_deg;

    endpoint_lat_max = std::max(endpoint_lat_max, lat_deg);
    endpoint_lat_min = std::min(endpoint_lat_min, lat_deg);
    endpoint_lon_max = std::max(endpoint_lon_max, lon_deg);
    endpoint_lon_min = std::min(endpoint_lon_min, lon_deg);
  }

  double lat_max =
      std::asin(std::clamp(info.face_z_max, -1.0, 1.0)) * rad_to_deg;
  double lat_min =
      std::asin(std::clamp(info.face_z_min, -1.0, 1.0)) * rad_to_deg;

  if (std::abs(lat_max - endpoint_lat_max) <= endpoint_lat_snap_tol_deg) {
    lat_max = endpoint_lat_max;
  }
  if (std::abs(lat_min - endpoint_lat_min) <= endpoint_lat_snap_tol_deg) {
    lat_min = endpoint_lat_min;
  }

  if (north_loc != Location::Outside) {
    bounds.lat_min = lat_min;
    bounds.lat_max = 90.0;
    bounds.lon_min = (north_loc == Location::Inside) ? 0.0 : endpoint_lon_min;
    bounds.lon_max = (north_loc == Location::Inside) ? 360.0 : endpoint_lon_max;
    bounds.longitude_wraps = (north_loc == Location::Inside);
    return bounds;
  }

  bounds.lat_min = -90.0;
  bounds.lat_max = lat_max;
  bounds.lon_min = (south_loc == Location::Inside) ? 0.0 : endpoint_lon_min;
  bounds.lon_max = (south_loc == Location::Inside) ? 360.0 : endpoint_lon_max;
  bounds.longitude_wraps = (south_loc == Location::Inside);
  return bounds;
}

inline LatLonBounds<double> generate_lat_lon_bounds(
    const double* const* face_vertices,
    std::size_t n,
    double polar_cap_lat_deg = default_polar_cap_lat_deg,
    double endpoint_lat_snap_tol_deg = default_endpoint_lat_snap_tol_deg) {
  const FaceLocationInfo<double> info =
      get_face_location_info(face_vertices, n, polar_cap_lat_deg);
  if (info.label == FaceLocationLabel::Local) {
    return generate_lat_lon_bounds_local(face_vertices, n,
                                         endpoint_lat_snap_tol_deg);
  }
  return generate_lat_lon_bounds_pole(face_vertices, n, polar_cap_lat_deg,
                                      endpoint_lat_snap_tol_deg);
}

namespace detail {

template <typename Container>
inline std::vector<const double*> to_lat_lon_ptr_vector(
    const Container& face_vertices) {
  std::vector<const double*> ptrs;
  ptrs.reserve(face_vertices.size());
  for (const auto& vertex : face_vertices) {
    ptrs.push_back(vertex.data());
  }
  return ptrs;
}

}  // namespace detail

template <std::size_t N>
inline FaceLocationInfo<double> get_face_location_info(
    const std::array<numeric::Vec3<double>, N>& face_vertices,
    double polar_cap_lat_deg = default_polar_cap_lat_deg) {
  const auto ptrs = detail::to_lat_lon_ptr_vector(face_vertices);
  return get_face_location_info(ptrs.data(), N, polar_cap_lat_deg);
}

inline FaceLocationInfo<double> get_face_location_info(
    const std::vector<numeric::Vec3<double>>& face_vertices,
    double polar_cap_lat_deg = default_polar_cap_lat_deg) {
  const auto ptrs = detail::to_lat_lon_ptr_vector(face_vertices);
  return get_face_location_info(ptrs.data(), face_vertices.size(),
                                polar_cap_lat_deg);
}

template <std::size_t N>
inline FaceLocationLabel get_face_location_label(
    const std::array<numeric::Vec3<double>, N>& face_vertices,
    double polar_cap_lat_deg = default_polar_cap_lat_deg) {
  return get_face_location_info(face_vertices, polar_cap_lat_deg).label;
}

inline FaceLocationLabel get_face_location_label(
    const std::vector<numeric::Vec3<double>>& face_vertices,
    double polar_cap_lat_deg = default_polar_cap_lat_deg) {
  return get_face_location_info(face_vertices, polar_cap_lat_deg).label;
}

template <typename T, std::size_t N>
inline FaceLocationInfo<T> get_face_location_info(
    const std::array<numeric::Vec3<T>, N>& face_vertices,
    T polar_cap_lat_deg = static_cast<T>(default_polar_cap_lat_deg)) {
  T face_z_max = -std::numeric_limits<T>::infinity();
  T face_z_min = std::numeric_limits<T>::infinity();

  for (std::size_t i = 0; i < N; ++i) {
    const numeric::Vec3<T>& x1 = face_vertices[i];
    const numeric::Vec3<T>& x2 = face_vertices[(i + 1) % N];

    const T z1 = x1[2];
    const T z2 = x2[2];
    const T d = x1[0] * x2[0] + x1[1] * x2[1] + x1[2] * x2[2];
    T a = (z1 * d - z2) / ((z1 + z2) * (d - T(1)));
    a = std::clamp(a, T(0), T(1));

    const T one_minus_a = T(1) - a;
    const T y0 = one_minus_a * x1[0] + a * x2[0];
    const T y1 = one_minus_a * x1[1] + a * x2[1];
    const T y2 = one_minus_a * x1[2] + a * x2[2];
    const T norm = std::sqrt(y0 * y0 + y1 * y1 + y2 * y2);
    const T za = y2 / norm;

    face_z_max = std::max(face_z_max, za);
    face_z_min = std::min(face_z_min, za);
  }

  const T polar_cap_z =
      std::sin(polar_cap_lat_deg * numeric::pi<T> / T(180));

  FaceLocationLabel label = FaceLocationLabel::Local;
  if (face_z_max >= polar_cap_z) {
    label = FaceLocationLabel::NorthPolarCapCandidate;
  } else if (face_z_min <= -polar_cap_z) {
    label = FaceLocationLabel::SouthPolarCapCandidate;
  }

  return {label, face_z_max, face_z_min};
}

template <typename T, std::size_t N>
inline LatLonBounds<T> generate_lat_lon_bounds_local(
    const std::array<numeric::Vec3<T>, N>& face_vertices,
    const FaceLocationInfo<T>& info,
    T endpoint_lat_snap_tol_deg =
        static_cast<T>(default_endpoint_lat_snap_tol_deg)) {
  LatLonBounds<T> bounds{};

  const T rad_to_deg = T(180) / numeric::pi<T>;

  T endpoint_lat_max = -std::numeric_limits<T>::infinity();
  T endpoint_lat_min = std::numeric_limits<T>::infinity();
  T endpoint_lon_max = -std::numeric_limits<T>::infinity();
  T endpoint_lon_min = std::numeric_limits<T>::infinity();

  for (std::size_t i = 0; i < N; ++i) {
    const numeric::Vec3<T>& x = face_vertices[i];
    const T z = std::clamp(x[2], T(-1), T(1));
    const T lat_deg = std::asin(z) * rad_to_deg;
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

  bounds.lat_max = lat_max;
  bounds.lat_min = lat_min;
  bounds.lon_max = endpoint_lon_max;
  bounds.lon_min = endpoint_lon_min;
  bounds.longitude_wraps = false;

  return bounds;
}

template <typename T, std::size_t N>
inline LatLonBounds<T> generate_lat_lon_bounds_local(
    const std::array<numeric::Vec3<T>, N>& face_vertices,
    T polar_cap_lat_deg = static_cast<T>(default_polar_cap_lat_deg),
    T endpoint_lat_snap_tol_deg =
        static_cast<T>(default_endpoint_lat_snap_tol_deg)) {
  const FaceLocationInfo<T> info =
      get_face_location_info(face_vertices, polar_cap_lat_deg);
  return generate_lat_lon_bounds_local(face_vertices, info,
                                       endpoint_lat_snap_tol_deg);
}

template <std::size_t N>
inline LatLonBounds<double> generate_lat_lon_bounds_pole(
    const std::array<numeric::Vec3<double>, N>& face_vertices,
    const FaceLocationInfo<double>& info,
    const std::int64_t* global_vertex_ids = nullptr,
    const std::int64_t* pole_query_id = nullptr,
    const std::int64_t* waypoint_id = nullptr,
    double endpoint_lat_snap_tol_deg = default_endpoint_lat_snap_tol_deg) {
  LatLonBounds<double> bounds{};

  const double rad_to_deg = 180.0 / numeric::pi<double>;

  constexpr std::int64_t kNorthPoleQueryId =
      std::numeric_limits<std::int64_t>::max() - 4;
  constexpr std::int64_t kSouthPoleQueryId =
      std::numeric_limits<std::int64_t>::max() - 3;
  constexpr std::int64_t kWaypointId =
      std::numeric_limits<std::int64_t>::max() - 2;

  constexpr double north_pole[3] = {0.0, 0.0, 1.0};
  constexpr double south_pole[3] = {0.0, 0.0, -1.0};

  std::array<const double*, N> face_ptrs{};
  for (std::size_t i = 0; i < N; ++i) {
    face_ptrs[i] = face_vertices[i].data();
  }

  std::array<std::int64_t, N> local_ids{};
  const std::int64_t* ids = global_vertex_ids;
  if (ids == nullptr) {
    for (std::size_t i = 0; i < N; ++i) {
      local_ids[i] = static_cast<std::int64_t>(i);
    }
    ids = local_ids.data();
  }

  const std::int64_t north_q_id =
      pole_query_id ? *pole_query_id : kNorthPoleQueryId;
  const std::int64_t south_q_id =
      pole_query_id ? *pole_query_id : kSouthPoleQueryId;

  const std::int64_t r_id = waypoint_id ? *waypoint_id : kWaypointId;
  (void)r_id;

  const Location north_loc = point_in_polygon_sphere(
      north_pole, north_q_id, face_ptrs.data(), ids, N);
  const Location south_loc = point_in_polygon_sphere(
      south_pole, south_q_id, face_ptrs.data(), ids, N);

  if (north_loc == Location::Outside && south_loc == Location::Outside) {
    return generate_lat_lon_bounds_local(face_vertices, info,
                                         endpoint_lat_snap_tol_deg);
  }

  double endpoint_lat_max = -std::numeric_limits<double>::infinity();
  double endpoint_lat_min = std::numeric_limits<double>::infinity();
  double endpoint_lon_max = -std::numeric_limits<double>::infinity();
  double endpoint_lon_min = std::numeric_limits<double>::infinity();

  for (std::size_t i = 0; i < N; ++i) {
    const auto& x = face_vertices[i];
    const double lat_deg = std::asin(std::clamp(x[2], -1.0, 1.0)) * rad_to_deg;
    const double lon_deg = std::atan2(x[1], x[0]) * rad_to_deg;

    endpoint_lat_max = std::max(endpoint_lat_max, lat_deg);
    endpoint_lat_min = std::min(endpoint_lat_min, lat_deg);
    endpoint_lon_max = std::max(endpoint_lon_max, lon_deg);
    endpoint_lon_min = std::min(endpoint_lon_min, lon_deg);
  }

  double lat_max = std::asin(std::clamp(info.face_z_max, -1.0, 1.0)) * rad_to_deg;
  double lat_min = std::asin(std::clamp(info.face_z_min, -1.0, 1.0)) * rad_to_deg;

  if (std::abs(lat_max - endpoint_lat_max) <= endpoint_lat_snap_tol_deg) {
    lat_max = endpoint_lat_max;
  }
  if (std::abs(lat_min - endpoint_lat_min) <= endpoint_lat_snap_tol_deg) {
    lat_min = endpoint_lat_min;
  }

  if (north_loc != Location::Outside) {
    bounds.lat_min = lat_min;
    bounds.lat_max = 90.0;
    bounds.lon_min = (north_loc == Location::Inside) ? 0.0 : endpoint_lon_min;
    bounds.lon_max = (north_loc == Location::Inside) ? 360.0 : endpoint_lon_max;
    bounds.longitude_wraps = (north_loc == Location::Inside);
    return bounds;
  }

  bounds.lat_min = -90.0;
  bounds.lat_max = lat_max;
  bounds.lon_min = (south_loc == Location::Inside) ? 0.0 : endpoint_lon_min;
  bounds.lon_max = (south_loc == Location::Inside) ? 360.0 : endpoint_lon_max;
  bounds.longitude_wraps = (south_loc == Location::Inside);
  return bounds;
}

template <std::size_t N>
inline LatLonBounds<double> generate_lat_lon_bounds_pole(
    const std::array<numeric::Vec3<double>, N>& face_vertices,
    double polar_cap_lat_deg = default_polar_cap_lat_deg,
    const std::int64_t* global_vertex_ids = nullptr,
    const std::int64_t* pole_query_id = nullptr,
    const std::int64_t* waypoint_id = nullptr,
    double endpoint_lat_snap_tol_deg = default_endpoint_lat_snap_tol_deg) {
  const FaceLocationInfo<double> info =
      get_face_location_info(face_vertices, polar_cap_lat_deg);
  return generate_lat_lon_bounds_pole(face_vertices, info, global_vertex_ids,
                                      pole_query_id, waypoint_id,
                                      endpoint_lat_snap_tol_deg);
}

template <std::size_t N>
inline LatLonBounds<double> generate_lat_lon_bounds(
    const std::array<numeric::Vec3<double>, N>& face_vertices,
    double polar_cap_lat_deg = default_polar_cap_lat_deg,
    const std::int64_t* global_vertex_ids = nullptr,
    const std::int64_t* pole_query_id = nullptr,
    const std::int64_t* waypoint_id = nullptr,
    double endpoint_lat_snap_tol_deg = default_endpoint_lat_snap_tol_deg) {
  const FaceLocationInfo<double> info =
      get_face_location_info(face_vertices, polar_cap_lat_deg);
  if (info.label == FaceLocationLabel::Local) {
    return generate_lat_lon_bounds_local(face_vertices, info,
                                         endpoint_lat_snap_tol_deg);
  }
  return generate_lat_lon_bounds_pole(face_vertices, info, global_vertex_ids,
                                      pole_query_id, waypoint_id,
                                      endpoint_lat_snap_tol_deg);
}

}  // namespace accusphgeom::algorithms
