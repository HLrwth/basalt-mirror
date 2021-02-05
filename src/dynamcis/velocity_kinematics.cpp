#include <basalt/dynamics/velocity_kinematics.h>

namespace basalt {

VelocityInverseKinematics::VelocityInverseKinematics()
  : start_t_ns(0),
    lin_var(0),
    ang_var(0){}
VelocityInverseKinematics::VelocityInverseKinematics(int64_t start_t_ns)
  : start_t_ns(start_t_ns),
    lin_var(0),
    ang_var(0){}

void VelocityInverseKinematics::estimateVelocity(const PoseVelState& state0,const PoseVelState& state1){

  //map to robot base frame
  //base frame: x forward, z up
  double dt = (state1.t_ns - state0.t_ns) * 1e-9;
  Sophus::SE3d T_w_o0 = state0.T_w_i * T_i_o;
  Sophus::SE3d T_w_o1 = state1.T_w_i * T_i_o;
  Sophus::SE3d T_o0_o1 = T_w_o0.inverse() * T_w_o1;

  //position
  Eigen::Vector2d pos_diff_2d = (T_o0_o1.translation()).segment<2>(0); //in o0 frame
  double trans_diff = pos_diff_2d.norm();

  //orientation
  Eigen::Vector3d rot_diff = (T_o0_o1.so3().inverse()).log();
  double yaw_diff = rot_diff(2); //rotate wrt. base frame. 

  //estimate
  double es_lin_vel; //forward
  double es_ang_vel;
  double es_gamma;
  
  //pure rotation 
  if(trans_diff<=std::numeric_limits<double>::min()){
    es_lin_vel = 0.0;
    es_ang_vel = yaw_diff / dt;
    es_gamma = 0.0;



    return;
  }
  
  //pure translation
  if(pos_diff_2d(1)<=std::numeric_limits<double>::min()){
    es_lin_vel = trans_diff/dt;
    es_ang_vel = 0.0; 
    es_gamma = 0.0;


    return;
  }

  //on arc
  Eigen::Vector2d center;
  center(0) = 0.0;
  center(1) = pos_diff_2d.squaredNorm()/(2.0 * pos_diff_2d(1));

  double radius = center(1);
  BASALT_ASSERT((pos_diff_2d-center).norm()<=std::numeric_limits<double>::min());
  BASALT_ASSERT(center(1)-pos_diff_2d(1)>0.0);
  double angle = atan2(pos_diff_2d(0),center(1)-pos_diff_2d(1));  
  double arc = radius * angle; 

  es_lin_vel = arc/dt;
  es_ang_vel = angle/dt;
  es_gamma = yaw_diff/dt - es_ang_vel;

  return;
}

Eigen::Vector3d VelocityInverseKinematics::residual(const PoseVelState& state0,const PoseVelState& state1,
                Mat3N* d_res_d_state0, Mat3N* d_res_d_state1){
                  

}





} // namespace basalt
