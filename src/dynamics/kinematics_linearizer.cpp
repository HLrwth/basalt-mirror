#include <basalt/vi_estimator/keypoint_vio.h>

namespace basalt{


void KeypointVioEstimator::linearizeKinematics(
      const AbsOrderMap& aom,Eigen::MatrixXd& abs_H, Eigen::VectorXd& abs_b,
      double& kin_error, double& extr_error, double& extr_dt_error,
      const Eigen::aligned_map<int64_t, PoseVelBiasExtrStateWithLin>& frame_states,
      const Eigen::aligned_map<int64_t, VelocityInverseKinematics>& kin_meas,
      const double& camera_cmd_timeoffset_weight, const double& camera_base_extr_weight){

  kin_error = 0; extr_error =0; extr_dt_error=0;

  for( const auto& kinematics : kin_meas){
    int64_t start_t = kinematics.second.get_start_t_ns();
    int64_t end_t = kinematics.second.get_end_t_ns();

    if (aom.abs_order_map.count(start_t) == 0 ||
    aom.abs_order_map.count(end_t) == 0)
      continue;

    const size_t start_idx = aom.abs_order_map.at(start_t).first;
    const size_t end_idx = aom.abs_order_map.at(end_t).first;

    auto start_state = frame_states.at(start_t);
    auto end_state = frame_states.at(end_t);

    VelocityInverseKinematics::Mat4N d_res_d_start, d_res_d_end;
    d_res_d_start.setZero(); d_res_d_end.setZero();
    VelocityInverseKinematics::Mat4N d_res_d_extr0, d_res_d_extr1;
    d_res_d_extr0.setZero(); d_res_d_extr1.setZero();
    Eigen::Vector4d d_res_d_toff0; Eigen::Vector4d d_res_d_toff1;
    d_res_d_toff0.setZero(); d_res_d_toff1.setZero();

    Eigen::Vector4d res = kinematics.second.residual(
        start_state.getStateLin(), end_state.getStateLin(), 
        &d_res_d_start, &d_res_d_end, &d_res_d_extr0, &d_res_d_extr1, 
        &d_res_d_toff0, &d_res_d_toff1);

    if (start_state.isLinearized() || end_state.isLinearized()) {
      res = kinematics.second.residual(
          start_state.getState(), end_state.getState());
    }

    Eigen::Matrix4d kinematics_cov_inv = kinematics.second.get_cov_inv();

    // error
    kin_error += 0.5 * res.transpose() * kinematics_cov_inv * res;

    // poses
    Eigen::Matrix<double,6,4> Jpose0weighted = d_res_d_start.transpose() * kinematics_cov_inv;
    Eigen::Matrix<double,6,4> Jpose1weighted = d_res_d_end.transpose() * kinematics_cov_inv;
    Eigen::Matrix<double,6,6> Jpose0pose1 = Jpose0weighted * d_res_d_end;
    abs_H.block<6 ,6>(start_idx, start_idx) += Jpose0weighted * d_res_d_start;
    abs_H.block<6 ,6>(start_idx, end_idx) += Jpose0pose1;
    abs_H.block<6 ,6>(end_idx, start_idx) += Jpose0pose1.transpose();
    abs_H.block<6 ,6>(end_idx, end_idx) += Jpose1weighted * d_res_d_end;

    abs_b.segment<6>(start_idx) += Jpose0weighted * res;
    abs_b.segment<6>(end_idx) += Jpose1weighted * res;

    // extrinsic
    Eigen::Matrix<double,6,4> Jextr0weighted = d_res_d_extr0.transpose() * kinematics_cov_inv;
    Eigen::Matrix<double,6,4> Jextr1weighted = d_res_d_extr1.transpose() * kinematics_cov_inv;
    Eigen::Matrix<double,6,6> Jextr0extr1 = Jextr0weighted * d_res_d_extr1;
    abs_H.block<6 ,6>(start_idx + 15, start_idx + 15) += Jextr0weighted * d_res_d_extr0;
    abs_H.block<6 ,6>(start_idx + 15, end_idx + 15) += Jextr0extr1;
    abs_H.block<6 ,6>(end_idx + 15, start_idx + 15) += Jextr0extr1.transpose();
    abs_H.block<6 ,6>(end_idx + 15, end_idx + 15) += Jextr1weighted * d_res_d_extr1;

    Eigen::Matrix<double,6,6> Jextr0pose0 = Jextr0weighted * d_res_d_start;
    abs_H.block<6 ,6>(start_idx + 15, start_idx) += Jextr0pose0;
    abs_H.block<6 ,6>(start_idx, start_idx + 15) += Jextr0pose0.transpose();

    Eigen::Matrix<double,6,6> Jextr0pose1 = Jextr0weighted * d_res_d_end;
    abs_H.block<6 ,6>(start_idx + 15, end_idx) += Jextr0pose1;
    abs_H.block<6 ,6>(end_idx, start_idx + 15) += Jextr0pose1.transpose();

    Eigen::Matrix<double,6,6> Jextr1pose0 = Jextr1weighted * d_res_d_start;
    abs_H.block<6 ,6>(end_idx + 15, start_idx) += Jextr1pose0;
    abs_H.block<6 ,6>(start_idx, end_idx + 15) += Jextr1pose0.transpose();

    Eigen::Matrix<double,6,6> Jextr1pose1 = Jextr1weighted * d_res_d_end;
    abs_H.block<6 ,6>(end_idx + 15, end_idx) += Jextr1pose1;
    abs_H.block<6 ,6>(end_idx, end_idx + 15) += Jextr1pose1.transpose();

    abs_b.segment<6>(start_idx + 15) += Jextr0weighted * res;
    abs_b.segment<6>(end_idx + 15) += Jextr1weighted * res;

    //time offset
    Eigen::Matrix<double,1,4> Jtoff0weighted = d_res_d_toff0.transpose() * kinematics_cov_inv;
    Eigen::Matrix<double,1,4> Jtoff1weighted = d_res_d_toff1.transpose() * kinematics_cov_inv;
    double Jtoff0toff1 = Jtoff0weighted.dot(d_res_d_toff1);
    abs_H(start_idx + 21, start_idx + 21) += Jtoff0weighted.dot(d_res_d_toff0);
    abs_H(start_idx + 21, end_idx + 21) += Jtoff0toff1;
    abs_H(end_idx + 21, start_idx + 21)+= Jtoff0toff1;
    abs_H(end_idx + 21, end_idx + 21)+= Jtoff1weighted.dot(d_res_d_toff1);

    Eigen::Matrix<double,1,6> Jtoff0pose0 = Jtoff0weighted * d_res_d_start;
    Eigen::Matrix<double,1,6> Jtoff0pose1 = Jtoff0weighted * d_res_d_end;
    abs_H.block<1 ,6>(start_idx + 21, start_idx) += Jtoff0pose0;
    abs_H.block<6, 1>(start_idx, start_idx + 21) += Jtoff0pose0.transpose();
    abs_H.block<1 ,6>(start_idx + 21, end_idx) += Jtoff0pose1;
    abs_H.block<6, 1>(end_idx, start_idx + 21) += Jtoff0pose1.transpose();

    Eigen::Matrix<double,1,6> Jtoff0extr0 = Jtoff0weighted * d_res_d_extr0;
    Eigen::Matrix<double,1,6> Jtoff0extr1 = Jtoff0weighted * d_res_d_extr1;
    abs_H.block<1 ,6>(start_idx + 21, start_idx + 15) += Jtoff0extr0;
    abs_H.block<6, 1>(start_idx + 15, start_idx + 21) += Jtoff0extr0.transpose();
    abs_H.block<1 ,6>(start_idx + 21, end_idx + 15) += Jtoff0extr1;
    abs_H.block<6, 1>(end_idx + 15, start_idx + 21) += Jtoff0extr1.transpose();

    Eigen::Matrix<double,1,6> Jtoff1pose0 = Jtoff1weighted * d_res_d_start;
    Eigen::Matrix<double,1,6> Jtoff1pose1 = Jtoff1weighted * d_res_d_end;
    abs_H.block<1 ,6>(end_idx + 21, start_idx) += Jtoff1pose0;
    abs_H.block<6 ,1>(start_idx, end_idx + 21) += Jtoff1pose0.transpose();
    abs_H.block<1 ,6>(end_idx + 21, end_idx) += Jtoff1pose1;
    abs_H.block<6, 1>(end_idx, end_idx + 21) += Jtoff1pose1.transpose();

    Eigen::Matrix<double,1,6> Jtoff1extr0 = Jtoff1weighted * d_res_d_extr0;
    Eigen::Matrix<double,1,6> Jtoff1extr1 = Jtoff1weighted * d_res_d_extr1;
    abs_H.block<1 ,6>(end_idx + 21, start_idx + 15) += Jtoff1extr0;
    abs_H.block<6, 1>(start_idx + 15, end_idx + 21) += Jtoff1extr0.transpose();
    abs_H.block<1 ,6>(end_idx + 21, end_idx + 15) += Jtoff1extr1;
    abs_H.block<6, 1>(end_idx + 15, end_idx + 21) += Jtoff1extr1.transpose();

    abs_b(start_idx + 21) += Jtoff0weighted * res;
    abs_b(end_idx + 21) += Jtoff1weighted * res;


    double dt = (end_t - start_t) * 1e-9;
    //difference between extrinsic
    {
      double camera_base_extr_weight_dt = camera_base_extr_weight * dt;
      Eigen::Matrix<double,6,1> res_extr;
      res_extr.head<3>() = end_state.getState().T_o_i.translation() - start_state.getState().T_o_i.translation();
      res_extr.tail<3>() = (start_state.getState().T_o_i.so3().inverse() * end_state.getState().T_o_i.so3()).log();

      Sophus::SO3d T_i0_o_Lin = start_state.getStateLin().T_o_i.so3().inverse();
      Eigen::Vector3d extrLin_rot_diff = (T_i0_o_Lin * end_state.getStateLin().T_o_i.so3()).log();
      Eigen::Matrix<double,3,3> J_left;
      Sophus::leftJacobianInvSO3(extrLin_rot_diff,J_left);

      Eigen::Matrix<double,6,6> d_diff_extr1;
      d_diff_extr1.setIdentity(); 
      d_diff_extr1.block<3, 3>(3, 3) = J_left * T_i0_o_Lin.matrix();

      Eigen::Matrix<double,6,6> J_diff_extr1weighted = d_diff_extr1.transpose() * camera_base_extr_weight_dt;
      Eigen::Matrix<double,6,6> J_diff_extr1extr1 = J_diff_extr1weighted * d_diff_extr1;
      abs_H.block<6 ,6>(start_idx + 15, start_idx + 15) += J_diff_extr1extr1;
      abs_H.block<6 ,6>(start_idx + 15, end_idx + 15) -= J_diff_extr1extr1;
      abs_H.block<6 ,6>(end_idx + 15, start_idx + 15) -= J_diff_extr1extr1;
      abs_H.block<6 ,6>(end_idx + 15, end_idx + 15) += J_diff_extr1extr1;

      abs_b.segment<6>(start_idx + 15) -= J_diff_extr1weighted * res_extr;
      abs_b.segment<6>(end_idx + 15) += J_diff_extr1weighted * res_extr;

      extr_error += 0.5 * res_extr.transpose() * camera_base_extr_weight_dt * res_extr;
    }

    //difference between timeoffset
    {
      double timeoffset_weight_dt = camera_cmd_timeoffset_weight * dt; 
      double res_timeoffset = start_state.getState().dt_extr - end_state.getState().dt_extr;

      abs_H(start_idx + 21 , start_idx + 21) += timeoffset_weight_dt;
      abs_H(end_idx + 21 , end_idx + 21) += timeoffset_weight_dt;

      abs_H(start_idx + 21, end_idx + 21) -= timeoffset_weight_dt;
      abs_H(end_idx + 21, start_idx + 21) -= timeoffset_weight_dt;

      abs_b(start_idx + 21) += timeoffset_weight_dt * res_timeoffset;
      abs_b(end_idx + 21) -= timeoffset_weight_dt * res_timeoffset;

      extr_dt_error += 0.5 * res_timeoffset * timeoffset_weight_dt * res_timeoffset;
    }

  }
}

void KeypointVioEstimator::computeKinematicsError(
      const AbsOrderMap& aom, double& kin_error, double& extr_error, double& extr_dt_error,
      const Eigen::aligned_map<int64_t, PoseVelBiasExtrStateWithLin>& frame_states,
      const Eigen::aligned_map<int64_t, VelocityInverseKinematics>& kin_meas,
      const double& camera_cmd_timeoffset_weight, const double& camera_base_extr_weight){
  kin_error=0; extr_error=0; extr_dt_error=0;

  for( const auto& kinematics : kin_meas){

    int64_t start_t = kinematics.second.get_start_t_ns();
    int64_t end_t = kinematics.second.get_end_t_ns();

    if (aom.abs_order_map.count(start_t) == 0 ||
    aom.abs_order_map.count(end_t) == 0)
      continue;

    auto start_state = frame_states.at(start_t);
    auto end_state = frame_states.at(end_t);

    const Eigen::Vector4d res = kinematics.second.residual(
          start_state.getState(), end_state.getState());
    
    kin_error += 0.5 * res.transpose() * kinematics.second.get_cov_inv() * res;

    double dt = (end_t - start_t) * 1e-9;
    {
      double camera_base_extr_weight_dt = camera_base_extr_weight / dt;
      Eigen::Matrix<double,6,1> res_extr;
      res_extr.head<3>() = end_state.getState().T_o_i.translation() - start_state.getState().T_o_i.translation();
      res_extr.tail<3>() = (start_state.getState().T_o_i.so3().inverse() * end_state.getState().T_o_i.so3()).log();

      extr_error += 0.5 * res_extr.transpose() * camera_base_extr_weight_dt * res_extr;
    }

    {
      double timeoffset_weight_dt = camera_cmd_timeoffset_weight / dt; 

      double res_timeoffset = start_state.getState().dt_extr - end_state.getState().dt_extr;
      extr_dt_error += 0.5 * res_timeoffset * timeoffset_weight_dt * res_timeoffset;
    }
  }
}



} //namesapce basalt