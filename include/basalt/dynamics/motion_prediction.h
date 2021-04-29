#pragma once

#include "basalt/dynamics/vel_command.h"
#include <basalt/imu/imu_types.h>

namespace basalt{

class MotionPrediction{

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  MotionPrediction(unsigned lag);

  void predict(const Eigen::aligned_map<int64_t, PoseVelBiasExtrStateWithLin>& frame_states,
               const Eigen::aligned_map<int64_t, VelocityInverseKinematics>& kin_meas); 





private:
  unsigned lag;
  VelCommand pre_vel_cmd;

  int64_t start_t_ns;
  int64_t end_t_ns;

  Sophus::SE3d T_w_i_pred;

  //tbb::concurrent_bounded_queue<basalt::PoseVelBiasState::Ptr> out_state_queue;

};









}//namespace