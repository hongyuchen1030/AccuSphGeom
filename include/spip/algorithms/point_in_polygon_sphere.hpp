#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace spip::pip {

enum class Location : std::uint8_t { Outside, Inside, OnVertex, OnEdge };

// No global IDs: use existing internal degeneracy rule.
Location point_in_polygon_sphere(const double* q,
                                 const double* const* poly,
                                 std::size_t n);

// With global IDs: use global-ID-based SoS perturbation rule.
Location point_in_polygon_sphere(const double* q,
                                 const double* const* poly,
                                 const std::int64_t* global_vertex_ids,
                                 std::size_t n);

Location point_in_polygon_sphere(const std::array<double, 3>& q,
                                 const std::vector<std::array<double, 3>>& poly);

Location point_in_polygon_sphere(const std::array<double, 3>& q,
                                 const std::vector<std::array<double, 3>>& poly,
                                 const std::vector<std::int64_t>& global_vertex_ids);

Location point_in_polygon_sphere(const std::vector<double>& q,
                                 const std::vector<std::vector<double>>& poly);

Location point_in_polygon_sphere(const std::vector<double>& q,
                                 const std::vector<std::vector<double>>& poly,
                                 const std::vector<std::int64_t>& global_vertex_ids);

}  // namespace spip::pip