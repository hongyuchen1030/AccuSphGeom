#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace spip::pip {

enum class Location : std::uint8_t { Outside, Inside, OnVertex, OnEdge };

Location point_in_polygon_sphere(const double* q,
                                 const double* const* poly,
                                 std::size_t n);

Location point_in_polygon_sphere(const std::array<double, 3>& q,
                                 const std::vector<std::array<double, 3>>& poly);

Location point_in_polygon_sphere(const std::vector<double>& q,
                                 const std::vector<std::vector<double>>& poly);

}  // namespace spip::pip