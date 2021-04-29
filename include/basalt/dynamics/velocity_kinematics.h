#pragma once

#include <basalt/dynamics/vel_command.h>
#include <basalt/imu/imu_types.h>

namespace basalt {

class VelocityInverseKinematics{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  using Ptr =  std::shared_ptr<VelocityInverseKinematics>;
  using Mat4N = Eigen::Matrix<double, 4, POSE_SIZE>;


  VelocityInverseKinematics();
  VelocityInverseKinematics(int64_t start_t_ns, const Eigen::Vector3d& start_gyro, 
                            int64_t end_t_ns , const Eigen::Vector3d& end_gyro, 
                            const VelCommand& vel_cmd);

  Eigen::Vector4d residual(const PoseVelBiasExtrState& state0,const PoseVelBiasExtrState& state1,
                Mat4N* d_res_d_state0 = nullptr, Mat4N* d_res_d_state1 = nullptr,
                Mat4N* d_res_d_extr = nullptr, Mat4N* d_res_d_extr1 = nullptr, 
                Eigen::Vector4d* d_res_d_toff0 = nullptr, Eigen::Vector4d* d_res_d_toff1 = nullptr) const;

  int64_t get_start_t_ns() const { return start_t_ns; }
  int64_t get_end_t_ns() const { return end_t_ns; }

  VelCommand get_vel_cmd() const { return vel_cmd; }

  inline const Eigen::Matrix4d get_cov_inv() const {
    if (!sqrt_cov_inv_computed) {
      compute_sqrt_cov_inv();
      sqrt_cov_inv_computed = true;
    }

    return sqrt_cov_inv.transpose() * sqrt_cov_inv;
  }

  /// @brief Square root inverse of the measurement covariance matrix
  inline const Eigen::Matrix4d& get_sqrt_cov_inv() const {
    if (!sqrt_cov_inv_computed) {
      compute_sqrt_cov_inv();
      sqrt_cov_inv_computed = true;
    }
    return sqrt_cov_inv;
  }

private:
  void compute_sqrt_cov_inv() const {
    sqrt_cov_inv.setIdentity();
    for (size_t i = 0; i < 4; i++) {
      if (cov(i,i) < std::numeric_limits<double>::min()) {
        sqrt_cov_inv(i,i) = 0.0;
      } else {
        sqrt_cov_inv(i,i) = 1.0 / sqrt(cov(i,i));
      }
    }
  }
  int64_t start_t_ns;
  Eigen::Vector3d start_gyro;
  int64_t end_t_ns;
  Eigen::Vector3d end_gyro;

  //cov 
  Eigen::Matrix4d cov;
  mutable Eigen::Matrix4d sqrt_cov_inv;
  mutable bool sqrt_cov_inv_computed;

  VelCommand vel_cmd;

};


}  // namespace basalt