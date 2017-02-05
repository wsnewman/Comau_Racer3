/* 
 * File:   comau_fk_ik.h
 * Author: wsn
 *
 * Created Feb 4, 2017
 * header for forward and inverse kinematics for Comau Racer3 robot
 * see: http://www.comau.com/Download/robot/racer3/Comau_Racer3_workingareas.pdf
 */

// NOTE: qvec for these fncs must be q in DH coords!!  use q_vec_comau(i) + DH_q_offsets(i)
// also, since DH z0 points DOWN, need to define a separate base frame w/
// z_base UP and x_base parallel to z0;
// i.e., Rotx(PI)= [1  0  0;
//                  0  -1 0;
//                  0  0 -1]

#ifndef COMAU_RACER3_FK_IK_H
#define	COMAU_RACER3_FK_IK_H
#include <ros/ros.h>
#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <string>
#include <math.h>

typedef Eigen::Matrix<double, 6, 1> Vectorq6x1;
//#include <boost/shared_ptr.hpp>
//#include <task_variables/library.h>

const int NJNTS=6;

// list DH params here
// structure and dimensions are remarkably similar to ABB IRB120
//assigned directions of positive rotation are odd
//assigning DH z axes consistent w/ Comau-defined positive rotations, with reference to fig
// http://www.comau.com/Download/robot/racer3/Comau_Racer3_workingareas.pdf
// z0 is DOWN
// z1 is OUT of page
// z2 is INTO page
// z3 is to RIGHT (forearm roll axis)
// z4 is OUT of page (wrist bend axis)
// z5 is to RIGHT (toolflange rotation axis)
// z6 is defined OUT of toolflange (antiparallel to z5)
const double DH_a1=0.050; //shoulder offset--different from IRB120
const double DH_a2=0.270;
const double DH_a3=0.050;
const double DH_a4=0.0;
const double DH_a5=0.0;
const double DH_a6=0.0;


const double DH_d1 = 0.365; //humerus length; <0, since z0 points DOWN
const double DH_d2 = 0.0;
const double DH_d3 = 0.0;
const double DH_d4 = 0.30594; //forearm also has joint axis pointing "IN", so negative d4
const double DH_d5 = 0.0;
const double DH_d6 = 0.080; //toolflange rotation also points "IN", so negative d6

//notes reference "mantis" position, "pos 9, calibration position"
//disregarding positive-rotation definitions (to be corrected later), CHOOSE DH z axes:
//choose to define DH z axes (w/rt calibration-position diagram) as:
//z0 UP, z1 OUT, z2 OUT, z3 to LEFT, z4 OUT, z5 to LEFT, z6 to LEFT

const double DH_alpha1 = -M_PI/2.0; //
  //DH x1 points to left...offset dir from z0 to z1
const double DH_alpha2 = 0; //shoulder and elbow z axes are defined parallel

// DH home angle has humerus FORWARD; i.e. q_dh3=0 (zero elbow ang) is same as calibration pose
// Comau home has humerus UP--> comau home is at q2_DH = -90deg 
const double DH_alpha3 = -M_PI/2.0; //elbow z OUT, forearm roll to LEFT--> -90 deg 

const double DH_alpha4 = M_PI/2.0; //forearm roll z3 to LEFT, wrist bend z4 OUT of page
                     //z3, z4 intersect, so x4 as ambiguous.  CHOOSE alpha = +90 deg,
                     //then, in calibration pose, x4 is UP, parallel to x3
                     //then comau home = q4_DH home (and forearm roll in home pose is consistent)

const double DH_alpha5 = -M_PI/2.0; // z5, toolflange spin, points LEFT 
                                    // z4 points out of page; z4,z5 intersect, so x5 is
                                    //ambiguous; CHOOSE alpha5 = -90 deg, so x5 points UP
                                    //in calibration pose, so zero wrist bend in calibration
                                    // pose has q5_DH = COMAU home = 0
const double DH_alpha6 = 0; //toolflange frame defined to point OUT from
                               // face of toolflange, parallel to z5 
 
//qvec_DH = sign(i)*qvec_comau(i) + DH_q_offsets(i);  see sign correction below
//i.e., at q_comau=0, what is q_DH? (q_DH = DH_q_offset)
const double DH_q_offset1 = 0.0;  //q1_DH is same as Comau home for base rotation
const double DH_q_offset2 = -M_PI/2.0; //at Comau home, q2_DH = -pi/2 for humerus pointing UP
const double DH_q_offset3 = -M_PI/2.0; //for elbow, Comau calibration pose=DH home, i.e. appears as 90deg elbow ang
const double DH_q_offset4 = 0.0; //at Comau home, q4_DH agrees w/ Comau angle
const double DH_q_offset5 = 0.0; //at Comau home, q5_DH agrees w/ Comau angle
const double DH_q_offset6 = 0.0; // home angle for toolflange x-axis can be defined consistent w/ Comau

const double deg2rad = M_PI/180.0;

//joint limits in Comau angles:
const double Comau_q_max1 = deg2rad*170;
const double Comau_q_max2 = deg2rad*135;
const double Comau_q_max3 = deg2rad*90;
const double Comau_q_max4 = deg2rad*179; //really, +/- 200deg; restrict to simplify IK
const double Comau_q_max5 = deg2rad*125;
//toolflange can rotate +/- 2700 deg!! simplify for IK
const double Comau_q_max6 = deg2rad*180; //deg2rad*400; //deliberately reduce this, to avoid excess solutions

const double Comau_q_min1 = -deg2rad*170;
const double Comau_q_min2 = -deg2rad*95;
const double Comau_q_min3 = -deg2rad*155; 
const double Comau_q_min4 = -deg2rad*179;
const double Comau_q_min5 = -deg2rad*125;
const double Comau_q_min6 = -deg2rad*180; //

const double DH_a_params[]={DH_a1,DH_a2,DH_a3,DH_a4,DH_a5,DH_a6};
const double DH_d_params[NJNTS] = {DH_d1, DH_d2, DH_d3, DH_d4, DH_d5, DH_d6};
const double DH_alpha_params[NJNTS] = {DH_alpha1, DH_alpha2, DH_alpha3, DH_alpha4, DH_alpha5, DH_alpha6};
const double DH_q_offsets[NJNTS] = {DH_q_offset1, DH_q_offset2, DH_q_offset3, DH_q_offset4, DH_q_offset5, DH_q_offset6};
//joint ranges, in Comau angles
const double q_lower_limits[NJNTS] = {Comau_q_min1, Comau_q_min2, Comau_q_min3, Comau_q_min4, Comau_q_min5, Comau_q_min6};
const double q_upper_limits[NJNTS] = {Comau_q_max1, Comau_q_max2, Comau_q_max3, Comau_q_max4, Comau_q_max5, Comau_q_max6};
//choose to define DH z axes (w/rt calibration-position diagram) as:
//z0 UP, z1 OUT, z2 OUT, z3 to LEFT, z4 OUT, z5 to LEFT, z6 to LEFT
// then qvec_DH = q_sign_corrections(i)*qvec_comau(i) + DH_q_offsets(i);
const double q_sign_corrections[NJNTS] = {-1,1,-1,-1, 1, -1};

class Comau_fwd_solver {
public:
    Comau_fwd_solver(); //constructor; //const hand_s& hs, const atlas_frame& base_frame, double rot_ang);
    //atlas_hand_fwd_solver(const hand_s& hs, const atlas_frame& base_frame);
    Eigen::Affine3d fwd_kin_solve(const Vectorq6x1& q_vec); // given vector of q angles, compute fwd kin
    Eigen::Affine3d fwd_kin_solve_q_comau(const Vectorq6x1& q_comau); //use q_vec per Comau defs for this fnc
    Eigen::Affine3d fwd_kin_solve_of_q_DH(const Vectorq6x1& q_vec); //use q_vec per DH defs for this fnc

    Eigen::Matrix4d get_wrist_frame();
    //Eigen::MatrixXd get_Jacobian(const Vectorq6x1& q_vec);
    Vectorq6x1 q_DH_of_q_comau(const Vectorq6x1& q_comau); //fncs to convert btwn Comau vars and DH vars
    Vectorq6x1 q_comau_of_q_DH(const Vectorq6x1& q_DH);
    bool fit_joints_to_range(Vectorq6x1 &qvec); //in COMAU joint coords!
    bool fit_joints_to_range_q_DH(Vectorq6x1 &qvec);   
    bool fit_q_to_range(double q_min, double q_max, double &q);      
private:  
    Eigen::Matrix4d fwd_kin_solve_(const Vectorq6x1& q_vec);
    Eigen::Matrix4d A_mats[6], A_mat_products[6], A_tool; // note: tool A must also handle diff DH vs URDF frame-7 xform
    Eigen::Matrix4d A_base;
    Eigen::MatrixXd Jacobian;
};

//IK solver--inherit fncs from fwd solver
class Comau_IK_solver: public Comau_fwd_solver {
public:
    Comau_IK_solver(); //constructor; 

    // return the number of valid solutions; actual vector of solutions will require an accessor function
    int ik_solve_DH_angs(Eigen::Affine3d const& desired_hand_pose); // given hand pose, compute vector DH-space jnt solns
    int ik_solve_comau_angs(Eigen::Affine3d const& desired_hand_pose,std::vector<Vectorq6x1> &q_comau_solns); 
    
    void get_DH_solns(std::vector<Vectorq6x1> &q_solns);

    //Eigen::MatrixXd get_Jacobian(const Vectorq6x1& q_vec);
private:
  
    std::vector<Vectorq6x1> q6dof_solns;
    std::vector<Vectorq6x1> q_solns_fit;
    Eigen::Matrix4d A_mats[6], A_mat_products[6], A_tool; // note: tool A must also handle diff DH vs URDF frame-7 xform
    double L_humerus;
    double L_forearm;
    double phi_elbow;
    //given desired flange pose, fill up solns for q1, q2, q3 based on wrist position
    bool compute_q123_solns(Eigen::Affine3d const& desired_hand_pose, std::vector<Vectorq6x1> &q_solns);
    //double solve_for_theta2(double q1,Eigen::Vector3d w_des);
    bool solve_for_theta2(Eigen::Vector3d w_wrt_1,double r_goal, double q2_solns[2]);    
    bool solve_for_theta3(Eigen::Vector3d w_wrt_1,double r_goal, double q3_solns[2]); 
    //alt approach: given q1 and q3, solve for corresponding q2; should have a unique soln, if reachable
    bool solve_for_shoulder_q(Eigen::Vector3d w_wrt_1, double q_elbow, double &q_shoulder);
    bool solve_spherical_wrist(Vectorq6x1 q_in,Eigen::Matrix3d R_des, std::vector<Vectorq6x1> &q_solns);    
    //Eigen::MatrixXd Jacobian;
};

#endif	/* IRB120_IK_H */

