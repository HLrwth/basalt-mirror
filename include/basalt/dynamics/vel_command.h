#pragma once

#include <basalt/utils/sophus_utils.hpp>
#include <memory>

namespace basalt {

///@brief Timestamped velocity commands
struct VelCommand
{
  using Ptr = std::shared_ptr<VelCommand>;

  int64_t t_ns;           ///< timestamp in nanoseconds
  double linear;
  double angular;

  double linear_diff;
  double angular_diff;

  /// @brief Default constructor with zero measurements.
  VelCommand() {
    t_ns = 0;
    linear = 0;
    angular = 0;
    linear_diff = 0;
    angular_diff = 0;
  }

  // VelCommand(int64_t t_ns, double linear, double angular):
  //   t_ns(t_ns),linear(linear),angular(angular){}  
};


}  // namespace basalt