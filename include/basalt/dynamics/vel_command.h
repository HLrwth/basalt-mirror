#pragma once

#include <basalt/utils/sophus_utils.hpp>
#include <memory>

namespace basalt {

///@brief Timestamped velocity commands
struct VelCommand
{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW 

  using Ptr = std::shared_ptr<VelCommand>;

  int64_t t_ns;           ///< timestamp in nanoseconds
  Eigen::Vector3d angular;  ///< Accelerometer measurement
  Eigen::Vector3d linear;   ///< Gyroscope measurement

  /// @brief Default constructor with zero measurements.
  VelCommand() {
    angular.setZero();
    linear.setZero();
  }

};


}  // namespace basalt