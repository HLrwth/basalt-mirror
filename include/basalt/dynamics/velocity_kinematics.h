#pragma once

#include <basalt/dynamics/vel_command.h>
#include <basalt/imu/imu_types.h>

namespace basalt {

class VelocityInverseKinematics{
EIGEN_MAKE_ALIGNED_OPERATOR_NEW
public:
  using Ptr =  std::shared_ptr<VelocityInverseKinematics>;
  using Mat3N = Eigen::Matrix<double, 3, POSE_VEL_SIZE>;


  VelocityInverseKinematics();
  VelocityInverseKinematics(int64_t start_t_ns);

  void estimateVelocity(const PoseVelState& state0,const PoseVelState& state1);
  Eigen::Vector3d residual(const PoseVelState& state0,const PoseVelState& state1,
                Mat3N* d_res_d_state0 = nullptr, Mat3N* d_res_d_state1 = nullptr);

private:
  int64_t start_t_ns;
  int64_t delta_t;

  //linear velocity factors
  double alpha_1;
  double alpha_2;

  //angular velocity factors
  double alpha_3;
  double alpha_4;

  //cov
  double lin_var;
  double ang_var;

  //extrinsic
  Sophus::SE3d T_i_o;
};


}  // namespace basalt