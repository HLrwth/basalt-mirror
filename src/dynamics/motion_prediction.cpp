#include "basalt/dynamics/motion_prediction.h"

namespace basalt{



void MotionPrediction::predict(const Eigen::aligned_map<int64_t, PoseVelBiasExtrStateWithLin>& frame_states,
                               const Eigen::aligned_map<int64_t, VelocityInverseKinematics>& kin_meas) {

  if(!kin_meas.empty()){
    auto kinematics_itr = kin_meas.begin();

    int64_t start_t = kinematics_itr->second.get_start_t_ns();
    int64_t end_t = kinematics_itr->second.get_end_t_ns();
    double dt = (end_t - start_t) * 1e-9; //ns to s

    auto start_state = frame_states.at(start_t);
    auto end_state = frame_states.at(end_t);

    double t0_offset = start_state.dt_extr;
    if(t0_offset <= 0.0){
      //linearized model
      Sophus::SE3d T_o0_o1;
      Eigen::Vector3d trans(dt * kinematics_itr->get_vel_cmd().linear, 0 ,0);
      T_o0_o1.translation() = trans;
      Eigen::Vector3d rot_vec(0, 0, dt * kinematics_itr->get_vel_cmd().angular_diff);
      T_o0_o1.so3() = Sophus::SO3d::exp(rot_vec);

      T_w_i_pred = T_w_i_pred * start_state.T_o_i.inverse() * T_o0_o1 * end_state.T_o_i();
    }
    else{
      Sophus::SE3d T_o0_o1;
      Eigen::Vector3d trans(t0_offset * pre_vel_cmd.linear + dt -  * kinematics_itr->get_vel_cmd().linear, 0 ,0);
      T_o0_o1.translation() = trans;
      Eigen::Vector3d rot_vec(0, 0, dt * kinematics_itr->get_vel_cmd().angular_diff);
      T_o0_o1.so3() = Sophus::SO3d::exp(rot_vec);
    }



    pre_vel_cmd = kinematics_itr->get_vel_cmd();
  }
}




} //namespace