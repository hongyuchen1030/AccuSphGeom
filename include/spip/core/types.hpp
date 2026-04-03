#pragma once

#include <Eigen/Dense>

namespace spip {

template <typename T>
using V3_T = Eigen::Matrix<T, 3, 1>;

template <typename T>
using Arc_T = Eigen::Matrix<T, 3, 2>;

}  // namespace spip