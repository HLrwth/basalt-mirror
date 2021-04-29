#include <basalt/dynamics/velocity_kinematics.h>

#include <sophus/se2.hpp>

namespace basalt {

VelocityInverseKinematics::VelocityInverseKinematics()
  : start_t_ns(0),end_t_ns(-1),
    sqrt_cov_inv_computed(false){
      cov.setIdentity();
      sqrt_cov_inv.setIdentity();
      start_gyro.setZero();
      end_gyro.setZero();
}

VelocityInverseKinematics::VelocityInverseKinematics(int64_t start_t_ns, const Eigen::Vector3d& start_gyro, 
                                                     int64_t end_t_ns , const Eigen::Vector3d& end_gyro, 
                                                     const VelCommand& vel_cmd)
  : start_t_ns(start_t_ns), start_gyro(start_gyro),
    end_t_ns(end_t_ns), end_gyro(end_gyro),
    sqrt_cov_inv_computed(false),
    vel_cmd(vel_cmd){

      //Initialize covariance
      cov.setZero();

      if(vel_cmd.linear_diff > 0){
        cov(0,0)+=  vel_cmd.linear_diff * 1e3;
        cov(2,2)+=  vel_cmd.linear_diff * 1e3;
      }

      if(vel_cmd.angular_diff > 0){
        cov(0,0)+= vel_cmd.angular_diff * 1e3;
        cov(1,1)+= vel_cmd.angular_diff * 1e3;
        cov(2,2)+= vel_cmd.angular_diff * 1e3;
      }

      cov(0,0) = std::max(cov(0,0) , 1e-2);
      cov(1,1) = std::max(cov(1,1) , 1e-2);
      cov(2,2) = std::max(cov(2,2) , 1e-2);
      cov(3,3) = 1e-2;
}


Eigen::Vector4d VelocityInverseKinematics::residual(
                const PoseVelBiasExtrState& state0,const PoseVelBiasExtrState& state1,
                Mat4N* d_res_d_state0, Mat4N* d_res_d_state1,
                Mat4N* d_res_d_extr0, Mat4N* d_res_d_extr1, 
                Eigen::Vector4d* d_res_d_toff0, Eigen::Vector4d* d_res_d_toff1) const{

  //If velocity cmd timestamp is later than frame, warp the frame after the cmd timestamp
  double t0_offset = state0.dt_extr;
  double t1_offset = state1.dt_extr;
  double frame_dt = (state1.t_ns - state0.t_ns) * 1e-9;
  double kin_dt = frame_dt + (t1_offset - t0_offset);

  //velocity interpolation
  Sophus::SE3d T_w_i0, T_w_i1;
  T_w_i0.so3() = state0.T_w_i.so3() * Sophus::SO3d::exp( start_gyro * t0_offset); 
  T_w_i0.translation() = state0.T_w_i.translation() + state0.vel_w_i * t0_offset;
  T_w_i1.so3() = state1.T_w_i.so3() * Sophus::SO3d::exp( end_gyro * t1_offset);
  T_w_i1.translation() = state1.T_w_i.translation() + state1.vel_w_i * t1_offset;

  Sophus::SO3d R_o0_w_sophus = state0.T_o_i.so3() * T_w_i0.so3().inverse();

  Sophus::SE3d T_o0_o1;
  T_o0_o1.so3() = R_o0_w_sophus * T_w_i1.so3() * state1.T_o_i.so3().inverse();

  Eigen::Vector3d trans_tmp = -(T_o0_o1.so3() * state1.T_o_i.translation());
  Eigen::Vector3d trans_tmp2 = R_o0_w_sophus * (T_w_i1.translation() - T_w_i0.translation()) + trans_tmp;
  T_o0_o1.translation() = trans_tmp2 + state0.T_o_i.translation();

  //position
  Eigen::Vector3d pos_diff = T_o0_o1.translation(); //in o0 frame

  //orientation
  Eigen::Vector3d rot_diff = (T_o0_o1.so3()).log(); //in o0 frame
  double yaw_diff = rot_diff(2); 

  //se2
  Sophus::SO2d yaw_so2(yaw_diff);
  double halftheta = 0.5 * yaw_diff;
  double halftheta_by_tan_of_halftheta;
  Eigen::Vector2d z = yaw_so2.unit_complex(); //cos(yaw_diff) + i * sin(yaw_diff)
  double real_minus_one = z.x() - 1.0;

  if (std::abs(real_minus_one) < Sophus::Constants<double>::epsilon()) {
    halftheta_by_tan_of_halftheta =
        1.0 - 1.0 / 12.0 * yaw_diff * yaw_diff;
  } else {
    halftheta_by_tan_of_halftheta = -(halftheta * z.y()) / (real_minus_one);
  }
  Eigen::Matrix2d V_inv;
  V_inv << halftheta_by_tan_of_halftheta, halftheta, -halftheta,
      halftheta_by_tan_of_halftheta;
  
  Eigen::Vector3d so2_log;
  so2_log.head<2>() = V_inv * pos_diff.head<2>();
  so2_log(2) = yaw_diff;
  Eigen::Vector3d so2_vel = so2_log / kin_dt;
  Eigen::Vector3d estimated(so2_vel(0),so2_vel(2),so2_vel(1));
  //res
  Eigen::Vector4d res;
  res(0) = estimated(0) - vel_cmd.linear;
  res(1) = estimated(1) - vel_cmd.angular;
  res(2) = estimated(2);
  res(3) = pos_diff(2); //z difference 0

  // std::cout<<"kin_dt: "<<kin_dt<<std::endl;
  // std::cout<<"Estimate linear vel: " << so2_vel(0)<<" "<<so2_vel(1)<<" vs cmd_vel linear: "<<vel_cmd.linear<<" and cmd linear diff: "<<vel_cmd.linear_diff<<"\n";
  // std::cout<<"Estimate angular vel: "<< so2_vel(2)<<" vs cmd_vel angualr: "<<vel_cmd.angular<<" and cmd angular diff: "<<vel_cmd.angular_diff<<"\n";
  // Eigen::Vector3d frame_ang = (state0.T_w_i.so3().inverse() * state1.T_w_i.so3()).log() / frame_dt;
  // std::cout<<"ang0: "<<start_gyro.transpose()<<" ang1: "<<end_gyro.transpose()<<" intrp: "<<frame_ang.transpose()<<"\n\n";
  // std::cout<<"extr t0: "<<state0.dt_extr<<" and extr t1: "<<state1.dt_extr<<"\n\n";

  if( d_res_d_state0 || d_res_d_state1 || d_res_d_extr0 || d_res_d_extr1 || d_res_d_toff0 ||d_res_d_toff1){
    //se2
    Eigen::Vector2d J_so2_pos_theta;
    {
    if (std::abs(real_minus_one) < 1e-10) {
      J_so2_pos_theta(0) = pos_diff(1) / 2.0 - yaw_diff * pos_diff(0) / 6.0;
      J_so2_pos_theta(1) = -pos_diff(0) / 2.0 - yaw_diff * pos_diff(1) / 6.0;
    }
    else{
      double sin_sq = z.y() * z.y();
      double real_minus_one_sq = real_minus_one * real_minus_one;
      J_so2_pos_theta(0) = ((-pos_diff(0) * z.y() - pos_diff(0) * z.x() * yaw_diff) / real_minus_one -
                            pos_diff(0) * sin_sq * yaw_diff / real_minus_one_sq + pos_diff(1)) / 2.0;
      J_so2_pos_theta(1) = ((-pos_diff(1) * z.y() - pos_diff(1) * z.x() * yaw_diff) / real_minus_one -
                            pos_diff(1) * sin_sq * yaw_diff / real_minus_one_sq - pos_diff(0)) / 2.0;
    }
    }
    
    //Jacobian from time offset compensation 
    Eigen::Matrix3d J_frame_ang_dt0, J_frame_ang_dt1;
    Sophus::leftJacobianSO3(start_gyro * t0_offset, J_frame_ang_dt0);
    Sophus::leftJacobianSO3(end_gyro * t1_offset, J_frame_ang_dt1);

    Eigen::Matrix3d r0 = state0.T_w_i.so3().matrix();
    J_frame_ang_dt0 = r0 * J_frame_ang_dt0;
    J_frame_ang_dt1 = state1.T_w_i.so3().matrix() * J_frame_ang_dt1;

    Eigen::Vector3d J_R0_dt0 = J_frame_ang_dt0 * start_gyro;
    Eigen::Vector3d J_R1_dt1 = J_frame_ang_dt1 * end_gyro;

    //Jacobian of res wrt. compensated states
    Eigen::Matrix3d R_o0_w = R_o0_w_sophus.matrix();
    Eigen::Matrix3d R_o0_o1 = T_o0_o1.so3().matrix();

    Eigen::Matrix3d tmp = Sophus::SO3d::hat(trans_tmp);  //B
    Eigen::Matrix3d tmp2 = Sophus::SO3d::hat(trans_tmp2);  //A

    Eigen::Matrix3d J;
    Sophus::leftJacobianInvSO3(rot_diff, J);
    Eigen::Matrix<double,1,3> J_yaw = J.row(2);

    Eigen::Matrix3d J_P_P0 = -R_o0_w; //PP0
    Eigen::Matrix3d J_P_R0 = tmp2 * R_o0_w; //PR0
    Eigen::Matrix<double,1,3> J_R_R0 = -J_yaw * R_o0_w; //RR0

    Eigen::Matrix3d J_P_R1 = -tmp * R_o0_w ; //PR1

    Eigen::Matrix<double,2,3> J_res02_P0 = V_inv * J_P_P0.topRows<2>();
    Eigen::Matrix<double,2,3> tmp3 = J_so2_pos_theta * J_R_R0;
    if(d_res_d_state0){
      d_res_d_state0->setZero();

      //J_P_p0 = J_P_P0;
      //derivative of lin wrt p0
      d_res_d_state0->block<1,3>(0,0) = J_res02_P0.row(0);
      //deivative of ang wrt p0 is 0
      d_res_d_state0->block<1,3>(2,0) = J_res02_P0.row(1);
      d_res_d_state0->block<1,3>(3,0) = J_P_P0.row(2);

      //J_P_r0 = J_P_R0;
      //derivative of lin wrt r0
      Eigen::Matrix<double,2,3> J_res02_R0 = V_inv * J_P_R0.topRows<2>() + tmp3;
      d_res_d_state0->block<1,3>(0,3) = J_res02_R0.row(0);
      //derivative of ang wrt r0
      d_res_d_state0->block<1,3>(1,3) = J_R_R0;
      //deivative of gamma wrt r0 is 0
      d_res_d_state0->block<1,3>(2,3) = J_res02_R0.row(1);
      d_res_d_state0->block<1,3>(3,3) = J_P_R0.row(2);

      d_res_d_state0->topRows<3>() /= kin_dt;
    }

    if(d_res_d_state1){
      d_res_d_state1->setZero();

      //Eigen::Matrix3d J_P_p1 = -J_P_P0;
      d_res_d_state1->block<1,3>(0,0) = -J_res02_P0.row(0);
      d_res_d_state1->block<1,3>(2,0) = -J_res02_P0.row(1);
      d_res_d_state1->block<1,3>(3,0) = -J_P_P0.row(2);

      //Eigen::Matrix3d J_P_r1 = J_P_R1;
      Eigen::Matrix<double,2,3> J_res02_R1 = V_inv * J_P_R1.topRows<2>() - tmp3;
      d_res_d_state1->block<1,3>(0,3) = J_res02_R1.row(0);
      d_res_d_state1->block<1,3>(1,3) = -J_R_R0;
      d_res_d_state1->block<1,3>(2,3) = J_res02_R1.row(1);
      d_res_d_state1->block<1,3>(3,3) = J_P_R1.row(2);

      d_res_d_state1->topRows<3>() /= kin_dt;
    }

    if(d_res_d_extr0){
      d_res_d_extr0->setZero();

      d_res_d_extr0->block<1,2>(0,0) = V_inv.row(0); //J_Px
      d_res_d_extr0->block<1,2>(2,0) = V_inv.row(1);
      (*d_res_d_extr0)(3,2) = 1.0;

      Eigen::Matrix<double,2,3> J_res02_extr_r0 = -V_inv * tmp2.topRows<2>() + J_so2_pos_theta * J_yaw;
      d_res_d_extr0->block<1,3>(0,3) = J_res02_extr_r0.row(0);
      d_res_d_extr0->block<1,3>(1,3) = J_yaw;
      d_res_d_extr0->block<1,3>(2,3) = J_res02_extr_r0.row(1);
      d_res_d_extr0->block<1,3>(3,3) = -tmp2.row(2);

      d_res_d_extr0->topRows<3>() /= kin_dt;
    }

    if(d_res_d_extr1){
      d_res_d_extr1->setZero();

      Eigen::Matrix<double,2,3> J_res02_extr_p1 = -V_inv * R_o0_o1.topRows<2>();
      d_res_d_extr1->block<1,3>(0,0) = J_res02_extr_p1.row(0);
      d_res_d_extr1->block<1,3>(2,0) = J_res02_extr_p1.row(1);
      d_res_d_extr1->block<1,3>(3,0) = -R_o0_o1.row(2);

      Eigen::Matrix<double,1,3> J_R_extr_rot = -J_yaw * R_o0_o1;
      Eigen::Matrix3d J_P_extr_rot = tmp * R_o0_o1;
      Eigen::Matrix<double,2,3> J_res02_extr_r1 = V_inv * J_P_extr_rot.topRows<2>() + J_so2_pos_theta * J_R_extr_rot;
      d_res_d_extr1->block<1,3>(0,3) = J_res02_extr_r1.row(0);
      d_res_d_extr1->block<1,3>(1,3) = J_R_extr_rot;
      d_res_d_extr1->block<1,3>(2,3) = J_res02_extr_r1.row(1);
      d_res_d_extr1->block<1,3>(3,3) = J_P_extr_rot.row(2);

      d_res_d_extr1->topRows<3>() /= kin_dt;
    }

    if(d_res_d_toff0){
      d_res_d_toff0->setZero();

      auto J_R_toff0 = J_R_R0 * J_R0_dt0;
      Eigen::Vector3d J_P_toff0 = J_P_P0 * state0.vel_w_i  + J_P_R0 * J_R0_dt0;
      Eigen::Vector2d J_res02_toff0 = V_inv * J_P_toff0.head<2>() + J_so2_pos_theta * J_R_toff0;
      (*d_res_d_toff0)(0) = J_res02_toff0(0);
      (*d_res_d_toff0)(1) = J_R_toff0;
      (*d_res_d_toff0)(2) = J_res02_toff0(1);
      (*d_res_d_toff0)(3) = J_P_toff0(2);

      d_res_d_toff0->topRows<3>() = (d_res_d_toff0->topRows<3>() + estimated) / kin_dt;
    }

    if(d_res_d_toff1){
      d_res_d_toff1->setZero();

      auto J_R_toff1 = -J_R_R0 * J_R1_dt1;
      Eigen::Vector3d J_P_toff1 = -J_P_P0 * state1.vel_w_i + J_P_R1 * J_R1_dt1;
      Eigen::Vector2d J_res02_toff1 = V_inv * J_P_toff1.head<2>() + J_so2_pos_theta * J_R_toff1;
      (*d_res_d_toff1)(0) = J_res02_toff1(0);
      (*d_res_d_toff1)(1) = J_R_toff1;
      (*d_res_d_toff1)(2) = J_res02_toff1(1);
      (*d_res_d_toff1)(3) = J_P_toff1(2);

      d_res_d_toff1->topRows<3>() = (d_res_d_toff1->topRows<3>() - estimated) / kin_dt;
    }
  }
  
  return res;
}


} // namespace basalt
