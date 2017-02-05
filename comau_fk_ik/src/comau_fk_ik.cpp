// comau_fk_ik implementation library; start w/ fwd kin

#include <comau_fk_ik/comau_fk_ik.h>

// function for use w/ both fwd and inv kin
// NOTE: q must be q in DH coords!!  use q_vec_comau(i) + DH_q_offsets(i)
// also, since DH z0 points DOWN, need to define a separate base frame w/
// z_base UP and x_base parallel to z0;

Eigen::Matrix4d compute_A_of_DH(double a, double d, double alpha, double q) {
    Eigen::Matrix4d A;
    Eigen::Matrix3d R;
    Eigen::Vector3d p;

    A = Eigen::Matrix4d::Identity();
    R = Eigen::Matrix3d::Identity();
    //ROS_INFO("compute_A_of_DH: a,d,alpha,q = %f, %f %f %f",a,d,alpha,q);

    double cq = cos(q);
    double sq = sin(q);
    double sa = sin(alpha);
    double ca = cos(alpha);
    R(0, 0) = cq;
    R(0, 1) = -sq*ca; //% - sin(q(i))*cos(alpha);
    R(0, 2) = sq*sa; //%sin(q(i))*sin(alpha);
    R(1, 0) = sq;
    R(1, 1) = cq*ca; //%cos(q(i))*cos(alpha);
    R(1, 2) = -cq*sa; //%	
    //%R(3,1)= 0; %already done by default
    R(2, 1) = sa;
    R(2, 2) = ca;
    p(0) = a * cq;
    p(1) = a * sq;
    p(2) = d;
    A.block<3, 3>(0, 0) = R;
    A.col(3).head(3) = p;
    return A;
}

//alt fnc, just takes q in DH joint space and index of frame, 0 to 5, and uses global DH vars

Eigen::Matrix4d compute_A_of_DH(int i, double q_dh) {
    Eigen::Matrix4d A;
    Eigen::Matrix3d R;
    Eigen::Vector3d p;
    double a = DH_a_params[i];
    double d = DH_d_params[i];
    double alpha = DH_alpha_params[i];
    double q = q_dh; //q_abb + DH_q_offsets[i];

    A = Eigen::Matrix4d::Identity();
    R = Eigen::Matrix3d::Identity();
    //ROS_INFO("compute_A_of_DH: a,d,alpha,q = %f, %f %f %f",a,d,alpha,q);

    double cq = cos(q);
    double sq = sin(q);
    double sa = sin(alpha);
    double ca = cos(alpha);
    R(0, 0) = cq;
    R(0, 1) = -sq*ca; //% - sin(q(i))*cos(alpha);
    R(0, 2) = sq*sa; //%sin(q(i))*sin(alpha);
    R(1, 0) = sq;
    R(1, 1) = cq*ca; //%cos(q(i))*cos(alpha);
    R(1, 2) = -cq*sa; //%	
    //%R(3,1)= 0; %already done by default
    R(2, 1) = sa;
    R(2, 2) = ca;
    p(0) = a * cq;
    p(1) = a * sq;
    p(2) = d;
    A.block<3, 3>(0, 0) = R;
    A.col(3).head(3) = p;
    return A;
}

//solve the eqn K = A*cos(q) + B*sin(q) for q; return "true" if at least one soln is valid

bool solve_K_eq_Acos_plus_Bsin(double K, double A, double B, std::vector<double> &q_solns) {
    double r, cphi, sphi, phi, gamma;
    double KTOL = 0.000001; //arbitrary tolerance
    r = sqrt(A * A + B * B);
    phi = atan2(B, A);
    q_solns.clear();
    if (fabs(K) > fabs(r)) {
        //ROS_WARN("K/r is too large for a cos/sin soln");
        return false; //if |K/r|>1, no solns
    }
    //could still have trouble w/ K=0...
    if (fabs(K) < KTOL) {
        //ROS_WARN("K is too small for A,B,K trig soln: user error? ");
        return false; //illegal use of this fnc
    }
    gamma = acos(K / r);
    double soln1 = phi + gamma;
    double soln2 = phi - gamma;
    q_solns.push_back(soln1);
    q_solns.push_back(soln2);
    return true;
}

Comau_fwd_solver::Comau_fwd_solver() { 

    ROS_INFO("fwd_solver constructor");
    // i.e., Rotx(PI)= [1  0  0;
//                  0  -1 0;
//                  0  0 -1]
    //A_base = compute_A_of_DH(0, 0, M_PI, 0);
    //std::cout<<"A_base = "<<std::endl;
    //std::cout<<A_base<<std::endl;
}

/*  IN CASE WANT JACOBIAN LATER...
Eigen::MatrixXd irb120_hand_fwd_solver::get_Jacobian(const Vectorq6x1& q_vec) {
    solve(q_vec);
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, 6);
    Eigen::MatrixXd J_ang = Eigen::MatrixXd::Zero(3, 6);
    Eigen::MatrixXd J_trans = Eigen::MatrixXd::Zero(3, 6);
    Eigen::MatrixXd zvecs = Eigen::MatrixXd::Zero(3, 6);
    Eigen::MatrixXd rvecs = Eigen::MatrixXd::Zero(3, 6);
    Eigen::Matrix4d Apalm = A_mat_products[7];
    Eigen::MatrixXd O_palm = Apalm.block<3, 1>(0, 3);
    Eigen::Matrix4d Ai;
    Eigen::MatrixXd zvec, rvec;
    Eigen::Vector3d t1, t2;
    for (int i = 0; i < 6; i++) {
        Ai = A_mat_products[i];
        zvec = Ai.block<3, 1>(0, 2); //%strip off z axis of each movable frame
        zvecs.block<3, 1>(0, i) = zvec; //%and store them
        rvec = O_palm - Ai.block<3, 1>(0, 3); //%vector from origin of i'th frame to palm 
        rvecs.block<3, 1>(0, i) = rvec;
        J_ang.block<3, 1>(0, i) = zvecs.block<3, 1>(0, i);

        t1 = zvecs.block<3, 1>(0, i);
        t2 = rvecs.block<3, 1>(0, i);
        J_trans.block<3, 1>(0, i) = t1.cross(t2);
    }

    J.block<3, 6>(0, 0) = J_trans;
    J.block<3, 6>(3, 0) = J_ang;
    if (is_lhand(hs_))return mirror_J_to_lhand(J);
    return J;
}

 */

Vectorq6x1 Comau_fwd_solver::q_DH_of_q_comau(const Vectorq6x1& q_comau) {
    Vectorq6x1 q_DH;
    for (int i=0;i<NJNTS;i++) {
        q_DH[i] = q_sign_corrections[i]*q_comau[i] + DH_q_offsets[i];
    }
    return q_DH;
}

Vectorq6x1 Comau_fwd_solver::q_comau_of_q_DH(const Vectorq6x1& q_DH) {
    Vectorq6x1 q_comau;
    for (int i=0;i<NJNTS;i++) {
        q_comau[i] = q_sign_corrections[i]*(q_DH[i] - DH_q_offsets[i]);
    }    
    return q_comau;
}

Eigen::Affine3d Comau_fwd_solver::fwd_kin_solve_q_comau(const Vectorq6x1& q_comau) {
    Vectorq6x1 q_DH;
    q_DH = q_DH_of_q_comau(q_comau);
    Eigen::Affine3d A;
    A = fwd_kin_solve_of_q_DH(q_DH);
}

Eigen::Affine3d Comau_fwd_solver::fwd_kin_solve_of_q_DH(const Vectorq6x1& q_vec) {
    Eigen::Matrix4d M;

    M = fwd_kin_solve_(q_vec);
    Eigen::Affine3d A(M);
    return A;
}

Eigen::Matrix4d Comau_fwd_solver::get_wrist_frame() {
    return A_mat_products[4];
}

//this fnc uses Comau joint angles
bool Comau_fwd_solver::fit_joints_to_range(Vectorq6x1 &q_comau) {
    bool fits=true;
    bool does_fit;
    double q;
    for (int i=0;i<NJNTS;i++) {
        q = q_comau[i];
        does_fit = fit_q_to_range(q_lower_limits[i],q_upper_limits[i],q);
        q_comau[i] = q;
        fits = fits&&does_fit;
    }
    if (fits)
        return true;
    else
        return false;
}

bool Comau_fwd_solver::fit_q_to_range(double q_min, double q_max, double &q) {
    while (q<q_min) {
        q+= 2.0*M_PI;
    }
    while (q>q_max) {
        q-= 2.0*M_PI;
    }    
    if (q<q_min)
        return false;
    else
        return true;
}

bool Comau_fwd_solver::fit_joints_to_range_q_DH(Vectorq6x1 &q_DH) {
    Vectorq6x1 q_comau;
    q_comau = q_comau_of_q_DH(q_DH); //convert to Comau joint angles
    fit_joints_to_range(q_comau); //test ranges w/rt Comau joint angle limits
    q_DH = q_DH_of_q_comau(q_comau); //convert back to DH angles; 

}

//return soln out to tool flange; would still need to account for tool transform for gripper
// ASSUMES q_vec IS EXPRESSED IN DH COORDS (requires translation from Comau joint angles)
Eigen::Matrix4d Comau_fwd_solver::fwd_kin_solve_(const Vectorq6x1& q_vec) {
    Eigen::Matrix4d A = Eigen::Matrix4d::Identity();
    //%compute A matrix from frame i to frame i-1:
    Eigen::Matrix4d A_i_iminusi;
    Eigen::Matrix3d R;
    Eigen::Vector3d p;
    for (int i = 0; i < NJNTS; i++) {
        //A_i_iminusi = compute_A_of_DH(DH_a_params[i],DH_d_params[i],DH_alpha_params[i], q_vec[i] + DH_q_offsets[i] );
        A_i_iminusi = compute_A_of_DH(i, q_vec[i]);
        A_mats[i] = A_i_iminusi;
        //std::cout<<"q_vec[i] = "<<q_vec[i]<<std::endl;
        //std::cout << "A_mats[" << i << "]:" << std::endl;
        //std::cout << A_mats[i] << std::endl;
    }

    A_mat_products[0] = A_mats[0];
        //std::cout << "A_mat_products[0]:" << std::endl;
        //std::cout << A_mat_products[0] << std::endl;      
    for (int i = 1; i < NJNTS; i++) {
        A_mat_products[i] = A_mat_products[i - 1] * A_mats[i];
        //std::cout << "A_mat_products[" << i << "]:" << std::endl;
        //std::cout << A_mat_products[i] << std::endl;        
    }

    //Eigen::Matrix4d A5;
    //A5 = A_mat_products[4];

    //std::cout<<" wrist position: "<<A5(0,3)<<", "<<A5(1,3)<<", "<<A5(2,3)<<std::endl;
    return A_mat_products[5]; //flange frame, expressed w/rt base frame w/ z_base UP
}

Comau_IK_solver::Comau_IK_solver() {
    //constructor: 
    L_humerus = DH_a_params[1];
    double L3 = DH_d_params[3];
    double A2 = DH_a_params[2];
    L_forearm = sqrt(A2 * A2 + L3 * L3);
    
    phi_elbow=acos((A2*A2+L_forearm*L_forearm-L3*L3)/(2.0*A2*L_forearm));
}

//int Comau_IK_solver::ik_solve(Eigen::Affine3d const& desired_hand_pose) {
int Comau_IK_solver::ik_solve_comau_angs(Eigen::Affine3d const& desired_hand_pose,std::vector<Vectorq6x1> &q_comau_solns) {
    ik_solve_DH_angs(desired_hand_pose);
    std::vector<Vectorq6x1> q_DH_solns;
    get_DH_solns(q_DH_solns);
    q_comau_solns.clear();
    Vectorq6x1 q_comau_soln;
    int n_DH_solns = q_DH_solns.size();
    for (int i=0;i<n_DH_solns;i++) {
        q_comau_soln = q_comau_of_q_DH(q_DH_solns[i]);
        if (fit_joints_to_range(q_comau_soln)) {
            q_comau_solns.push_back(q_comau_soln);
        }
    }
    //Vectorq6x1 Comau_fwd_solver::q_comau_of_q_DH(const Vectorq6x1& q_DH) {
    int n_comau_solns = q_comau_solns.size();
    return n_comau_solns;
}


int Comau_IK_solver::ik_solve_DH_angs(Eigen::Affine3d const& desired_hand_pose) {    
    q6dof_solns.clear();
    bool reachable = compute_q123_solns(desired_hand_pose, q6dof_solns);
    if (!reachable) {
        return 0;
    }
    reachable = false;
    //is at least one solution within joint range limits?
    q_solns_fit.clear();
    Vectorq6x1 q_soln;
    Eigen::Matrix3d R_des;
    R_des = desired_hand_pose.linear();
    int nsolns = q6dof_solns.size();
    bool fits;

    std::vector<Vectorq6x1> q_wrist_solns;
    for (int i=0;i<nsolns;i++) {
        q_soln = q6dof_solns[i];
        //ROS_WARN("q1,q2,q3 before fit = %f, %f, %f",q_soln[0],q_soln[1],q_soln[2]);
        fits = fit_joints_to_range_q_DH(q_soln); // force q_soln in to periodic range, if possible, and return if possible
        if (fits) { // if here, then have a valid 3dof soln; try to get wrist solutions
            // get wrist solutions; expect 2, though not checked for joint limits
            //ROS_WARN("get wrist solns for q123[%d]",i);
            //ROS_WARN("q1,q2,q3 = %f, %f, %f",q_soln[0],q_soln[1],q_soln[2]);
            solve_spherical_wrist(q_soln,R_des, q_wrist_solns);  
            int n_wrist_solns = q_wrist_solns.size();
            for (int iwrist=0;iwrist<n_wrist_solns;iwrist++) {
                q_soln = q_wrist_solns[iwrist];
                if (fit_joints_to_range_q_DH(q_soln)) {
                  q_solns_fit.push_back(q_soln);
                  reachable = true; // note that we have at least one reachable solution
                }
            }

        }
    }
    if (!reachable) {
        return 0;
    }
    
    //if here, have reachable solutions
    nsolns = q_solns_fit.size();
    return nsolns;
}


//accessor function to get all solutions

void Comau_IK_solver::get_DH_solns(std::vector<Vectorq6x1> &q_solns) {
    q_solns = q_solns_fit; //q6dof_solns;
}

//given desired flange pose, fill up solns for q1, q2, q3 based on wrist position
// these angles are NOT checked for joint limits
bool Comau_IK_solver::compute_q123_solns(Eigen::Affine3d const& desired_hand_pose, std::vector<Vectorq6x1> &q_solns) {
    double L6 = DH_d_params[5];
    double r_goal;
    bool reachable=false;
    bool reachable_q1a = true;
    bool reachable_q1b = true;
    Eigen::Vector3d p_des = desired_hand_pose.translation();
    Eigen::Matrix3d R_des = desired_hand_pose.linear();
    Eigen::Vector3d z_des = R_des.col(2); // direction of desired z-vector
    Eigen::Vector3d w_des = p_des - L6*z_des; // desired wrist position w/rt frame0
    //std::cout<<"IK inferred des wrist point: "<<w_des.transpose()<<std::endl;
    q_solns.clear();
    Vectorq6x1 q_soln;
    
    double q1a = atan2(w_des(1), w_des(0));
    double q1b = q1a + M_PI; // given q1, q1+pi is also a soln
    Eigen::Matrix4d A1a, A1b;
    //compute_A_of_DH(double a,double d,double q, double alpha); use q_vec(i) + DH_q_offsets(i)
    A1a = compute_A_of_DH(0, q1a); //compute_A_of_DH(DH_a_params[0],DH_d_params[0],DH_alpha_params[0], q1a+ DH_q_offsets[0] );
    A1b = compute_A_of_DH(0, q1b); //compute_A_of_DH(DH_a_params[0],DH_d_params[0],DH_alpha_params[0], q1b+ DH_q_offsets[0] );    
    double q2a_solns[2];
    double q2b_solns[2];

    Eigen::Matrix4d A10;
    A10 = compute_A_of_DH(0, q1a);
    Eigen::Matrix3d R10;
    Eigen::Vector3d p1_wrt_0, w_wrt_1a, w_wrt_1b;
    R10 = A10.block<3, 3>(0, 0);
    //std::cout << "compute_q123, R10: " << std::endl;
    //std::cout << R10 << std::endl;
    p1_wrt_0 = A10.col(3).head(3); ///origin of frame1 w/rt frame0
    //compute A10_inv * w_wrt_0, = R'*w_wrt_0 -R'*p1_wrt_0
    w_wrt_1a = R10.transpose() * w_des - R10.transpose() * p1_wrt_0; //desired wrist pos w/rt frame1
    //std::cout << "w_wrt_1a = " << w_wrt_1a.transpose() << std::endl;
    r_goal = sqrt(w_wrt_1a[0] * w_wrt_1a[0] + w_wrt_1a[1] * w_wrt_1a[1]);

    //ROS_INFO("r_goal = %f", r_goal);
    //is the desired wrist position reachable? if not, return false
    // does not yet consider joint limits
    if (r_goal >= L_humerus + L_forearm) {
        //ROS_INFO("goal is too far away!");
        reachable_q1a= false;
    }
    // can also have problems if wrist goal is too close to shoulder
    if (r_goal <= fabs(L_humerus - L_forearm)) {
        //ROS_INFO("goal is too close!");
        reachable_q1a= false;
    }

    //reachable = solve_for_theta2(w_wrt_1a, r_goal, q2a_solns);
    double q3a_solns[2];
    if(reachable_q1a) {
      reachable_q1a = solve_for_theta3(w_wrt_1a, r_goal, q3a_solns);
    }

    if (!reachable_q1a) {
        //ROS_WARN(" wrist point is out of reach for q1a");
        return false;
    } 
    else {
        //ROS_INFO("q3a solns: %f, %f", q3a_solns[0], q3a_solns[1]);
    }    
    
    //ROS_INFO("try these solns as inputs for q_shoulder: ");
    double q_shoulder;
    double q_shoulder_solnsa[2];
    if (reachable_q1a) {
      bool q_shoulder_ok = solve_for_shoulder_q(w_wrt_1a, q3a_solns[0], q_shoulder);
      if (q_shoulder_ok) {
                //ROS_INFO("shoulder soln from first elbow ang: %f",q_shoulder);
              q_soln[0] = q1a;
              q_soln[1] = q_shoulder;
              q_soln[2] = q3a_solns[0];
              q_solns.push_back(q_soln);
              reachable=true; // have at least one soln
      }
      //q_shoulder_solnsa[0] = q_shoulder;

      q_shoulder_ok = solve_for_shoulder_q(w_wrt_1a, q3a_solns[1], q_shoulder);
      if (q_shoulder_ok) {
        //ROS_INFO("shoulder soln from second elbow ang: %f",q_shoulder);
        q_shoulder_solnsa[1] = q_shoulder;
              q_soln[0] = q1a;
              q_soln[1] = q_shoulder;
              q_soln[2] = q3a_solns[1];
              q_solns.push_back(q_soln);
              reachable=true; // have at least one soln        
      }
    }
    // now, get w_wrt_1b:
    //ROS_WARN("repeat for q1b");
    A10 = compute_A_of_DH(0, q1b);
    R10 = A10.block<3, 3>(0, 0);
    p1_wrt_0 = A10.col(3).head(3); //origin of frame1 w/rt frame0; should be the same as above
    //compute A10_inv * w_wrt_0, = R'*w_wrt_0 -R'*p1_wrt_0
    w_wrt_1b = R10.transpose() * w_des - R10.transpose() * p1_wrt_0; //desired wrist pos w/rt frame1
    //std::cout << "w_wrt_1b = " << w_wrt_1b.transpose() << std::endl;
    //r_goal is not same as above, due to a1 offset
    r_goal = sqrt(w_wrt_1b[0] * w_wrt_1b[0] + w_wrt_1b[1] * w_wrt_1b[1]);
    if (r_goal >= L_humerus + L_forearm) {
        //ROS_INFO("goal is too far away!");
        reachable_q1b= false;
    }
    // can also have problems if wrist goal is too close to shoulder
    if (r_goal <= fabs(L_humerus - L_forearm)) {
        //ROS_INFO("goal is too close!");
        reachable_q1b= false;
    }    
    //ROS_INFO("r_goal b: %f", r_goal);

    //reachable = solve_for_theta2(w_wrt_1b, r_goal, q2b_solns);
    //ROS_INFO("q2b solns: %f, %f", q2b_solns[0], q2b_solns[1]);
    
    double q3b_solns[2];
    if(reachable_q1b) {
      reachable_q1b = solve_for_theta3(w_wrt_1b, r_goal, q3b_solns);
    } 
    if (reachable_q1b) {
      bool q_shoulder_ok = solve_for_shoulder_q(w_wrt_1b, q3b_solns[0], q_shoulder);
      if (q_shoulder_ok) {
               // ROS_INFO("shoulder soln from first elbow ang: %f",q_shoulder);
              q_soln[0] = q1b;
              q_soln[1] = q_shoulder;
              q_soln[2] = q3b_solns[0];
              q_solns.push_back(q_soln);
              reachable=true; // have at least one soln
      }    
    
      q_shoulder_ok = solve_for_shoulder_q(w_wrt_1b, q3b_solns[1], q_shoulder);
      if (q_shoulder_ok) {
        //ROS_INFO("shoulder soln from second elbow ang: %f",q_shoulder);
        //q_shoulder_solnsb[1] = q_shoulder;
              q_soln[0] = q1b;
              q_soln[1] = q_shoulder;
              q_soln[2] = q3b_solns[1];
              q_solns.push_back(q_soln);
              reachable=true; // have at least one soln        
      }
    }
 
    
    return reachable;

}

//reachable = solve_for_theta2(q1a,w_wrt_1,r_goal,q2a_solns);
//solve for shoulder angle
bool Comau_IK_solver::solve_for_theta2(Eigen::Vector3d w_wrt_1, double r_goal, double q2_solns[2]) {
    //double L_humerus = DH_a_params[1];
    //double L3 = DH_d_params[3];
    //double A2 = DH_a_params[2];
    //double L_forearm = sqrt(A2*A2 + L3*L3);


    double beta = atan2(-w_wrt_1[1], w_wrt_1[0]);
    double acos_arg = (L_humerus * L_humerus + r_goal * r_goal - L_forearm * L_forearm) / (2.0 * r_goal * L_humerus);
    if (fabs(acos_arg > 1.0)) {
        ROS_WARN("hey!  logic err acos_arg = %f", acos_arg);
        return false;
    }
    double gamma = acos(acos_arg);
    q2_solns[0] = M_PI / 2.0 - beta - gamma; //%elbow up
    q2_solns[1] = M_PI / 2.0 - beta + gamma; //%elbow down 
    return true;
    //return 0;
}



//find elbow angle solutions
bool Comau_IK_solver::solve_for_theta3(Eigen::Vector3d w_wrt_1, double r_goal, double q3_solns[2]) {
     //phi_elbow;
    double acos_arg = (L_humerus*L_humerus + L_forearm*L_forearm - r_goal*r_goal)/(2.0*L_humerus*L_forearm);
    if (fabs(acos_arg>1.0)) {
        ROS_WARN("solve_for_theta3 logic err!  acos_arg = %f",acos_arg);
        return false;
    }
    double eta = acos(acos_arg);     
    
    q3_solns[0]=  M_PI -phi_elbow -eta;    //q3(1) = pi - phi- eta; 
    q3_solns[1]=  M_PI -phi_elbow +eta;  //q3(2) = pi - phi + eta; */
    return true;
}

//alt approach: given q1 and q3, solve for corresponding q2; should have a unique soln, if reachable
// use q1 to provide wrist point w/rt frame 1 == w_wrt_1
// provide q_elbow; solve for q_shoulder
//this fnc does not test joint limits
// ASSUMES shoulder angle in DH coords
bool Comau_IK_solver::solve_for_shoulder_q(Eigen::Vector3d w_wrt_1, double q_elbow, double &q_shoulder) {
    //double L_humerus = DH_a_params[1];
    //double L3 = DH_d_params[3];
    //double A2 = DH_a_params[2];
    //double L_forearm = sqrt(A2*A2 + L3*L3);
    //std::cout<<"trying to reach w_wrt_1 = "<<w_wrt_1.transpose()<<" at q_elbow = "<<q_elbow<<std::endl;
    //compute A32(q_elbow):
    Eigen::Matrix4d A43,A32,A21,A41,A42;
    //we know [wrist_pt;1] = A10(q1)*A21(q2)*A32(q3)*A43(q4)*[0;0;0;1]
    //note that A43(q4)*[0;0;0;1] = [0;0;d4;1], and does not depend on q4
    //also, given q1, compute A10 and premultiply by inverse:
    // A10_inv * [wrist_pt;1] = A21(q2) * A32(q3)*[0;0;d4;1]
    // given q3, can compute vector A32(q3)*[0;0;d4;1], call this vec324
    // note that vec324_z = 0; 
    // also, LHS is same as wrist point w/rt frame 1, which is provided as an argument
    A43 = compute_A_of_DH(3, 0); //transform frame 4 wrt frame 3 at forearm roll angle =0
    A32 = compute_A_of_DH(2, q_elbow); //transform A-matrix of frame 3 w/rt frame 2
    A42 = A32*A43;
    //A42 = A32*A43; extract last row of A42, which is A32*[0;0;d4;1]
    //here is the result, expressed as a 3-vector
    // geometrically, this is the point of intersection from J3 to J4 axes, as expressed in DH frame2
    Eigen::Vector3d vec324; // = DH_d_params[3]*A32.col(2).head(3)+A32.col(3).head(3);
    //vec324 = DH_d_params[3]*A32.col(2).head(3)+A32.col(3).head(3);
    vec324= A42.col(3).head(3);
    // use first row of w_wrt_1[0] = A21(0,0)*vec324(0) + A21(0,1)*vec324(1) + A21(0,2)*0 + A21(0,3)*1
    //i.e., w_x = cos(q2)*vec324_x -sin(q2)*vec324_y + a2*cos(q2)
    //or, w_x = (vec324_x + a2)*cos(q2) + (-vec324_y)*sin(q2)
    // solve K = A*cos(q) + B*sin(q)
    //use: bool solve_K_eq_Acos_plus_Bsin(double K, double A, double B, std::vector<double> &q_solns)
    std::vector<double> q_shoulder_solns;
    double K = w_wrt_1[0];
    double A = vec324[0] + DH_a_params[1];
    double B = -vec324[1];
    bool soln_ok = solve_K_eq_Acos_plus_Bsin(K,  A,  B, q_shoulder_solns);
    if (!soln_ok) return false;
    //test the two solutions to find which one applies:
    q_shoulder = q_shoulder_solns[0]; //default assumption
    A21 = compute_A_of_DH(1, q_shoulder_solns[0]);
    A41 = A21*A32*A43;
    Eigen::Vector3d w_test;
    w_test= A41.col(3).head(3); //extract wrist position w/rt frame 1 based on first q2 soln
    //std::cout<<"first wrist soln: "<<w_test.transpose()<<std::endl;
    double w_err = (w_test-w_wrt_1).norm();
    //ROS_INFO("First soln err = %f",w_err);
    A21 = compute_A_of_DH(1, q_shoulder_solns[1]);
    A41 = A21*A32*A43;   
    w_test= A41.col(3).head(3);
    //std::cout<<"second wrist soln: "<<w_test.transpose()<<std::endl;
    double w_err2 = (w_test-w_wrt_1).norm();
    //ROS_INFO("Second soln err = %f",w_err2);
    if (w_err2<w_err) q_shoulder = q_shoulder_solns[1];
    
    return true;
    //return 0;
}




//find the + forearm-rotation orientation solution, q4, q5, q6--not concerned with joint limits;
// note: if q5 is near zero, then at a wrist singularity; 
// inf solutions of q4+D, q6-D
// use q1, q2, q3 from q_in; copy these values to q_solns, and tack on the two solutions q4, q5, q6
//USES DH angles, not Comau angles
bool Comau_IK_solver::solve_spherical_wrist(Vectorq6x1 q_in,Eigen::Matrix3d R_des, std::vector<Vectorq6x1> &q_solns) {
    bool is_singular = false;
    Eigen::Matrix4d A01,A12,A23,A03,A34,A04,A45,A05;
    A01 = compute_A_of_DH(0, q_in[0]);
    A12 = compute_A_of_DH(1, q_in[1]);
    A23 = compute_A_of_DH(2, q_in[2]);
    A03 = A01*A12*A23;   
    Eigen::Vector3d n3,t3,b3; //axes of frame3
    Eigen::Vector3d n4,t4,b4; // axes of frame4
    Eigen::Vector3d n5,t5; // axes of frame5; b5 is same as b_des = b6
    Eigen::Vector3d n_des,b_des; // desired x-axis and z-axis of flange frame
    n3 = A03.col(0).head(3);
    t3 = A03.col(1).head(3);    
    b3 = A03.col(2).head(3); 
    //ROS_INFO("q1,q2,q3 = %f, %f, %f",q_in[0],q_in[1],q_in[2]);

    //std::cout<<"n3 = "<<n3.transpose()<<std::endl;
    //std::cout<<"t3 = "<<t3.transpose()<<std::endl;
    //std::cout<<"b3 = "<<b3.transpose()<<std::endl;
    
    b_des = R_des.col(2);
    n_des = R_des.col(0);
    b4 = b3.cross(b_des);
    double q4,q5,q6;
    Vectorq6x1 q_soln;
      if (b4.norm() <= 0.000001) {
                q4=0;
                is_singular = true;
      }
      else {
            double cq4= b4.dot(-t3);
            double sq4= b4.dot(n3); 
            q4= atan2(sq4, cq4);
        }
    // choose the positive forearm-rotation solution:
    if (q4>M_PI) {
        q4-= 2*M_PI;
    }    
    if (q4<0.0) {
        q4+= M_PI;
    }
    double q4b = q4 -M_PI;
    // THESE OPTIONS LOOK GOOD FOR q4
    //std::cout<<"forearm rotation options: "<<q4<<", "<<q4b<<std::endl;
    
    // use the + q4 soln to find q5, q6
    A34 = compute_A_of_DH(3, q4);
    A04 = A03*A34;
    //std::cout<<"A04 = "<<std::endl;
    //std::cout<<A04<<std::endl;
    n4 = A04.col(0).head(3);
    t4 = A04.col(1).head(3); 
    double cq5 = b_des.dot(t4);
    double sq5 = b_des.dot(-n4);
    q5 = atan2(sq5,cq5);
    //std::cout<<"bdes = "<<b_des.transpose()<<std::endl;
    //std::cout<<"at q4 = "<<q4<<std::endl;
    //std::cout<<"n4 = "<<n4.transpose()<<std::endl;
    //std::cout<<"wrist bend = "<<q5<<std::endl;

    //solve for q6
    A45 = compute_A_of_DH(4, q5);
    A05 = A04*A45;
    n5 = A05.col(0).head(3);
    t5 = A05.col(1).head(3);   
    //std::cout<<"ndes = "<<n_des.transpose()<<std::endl;
    //std::cout<<"n5 = "<<n5.transpose()<<std::endl; 
    //std::cout<<"t5 = "<<t5.transpose()<<std::endl;     
    double cq6=n_des.dot(n5);
    double sq6=n_des.dot(t5);
    q6 =atan2(sq6, cq6);
    //ROS_INFO("q4,q5,q6 = %f, %f, %f",q4,q5,q6);
    q_soln = q_in;
    q_soln[3] = q4;
    q_soln[4] = q5;
    q_soln[5] = q6;
    q_solns.clear();
    q_solns.push_back(q_soln);
    //2nd wrist soln: 
    q_soln[3] = q4b;
    q_soln[4] *= -1.0; // flip wrist opposite direction
    q_soln[5] = q6+M_PI; // fix the periodicity later; 
       // ROS_INFO("alt q4,q5,q6 = %f, %f, %f",q_soln[3],q_soln[4],q_soln[5]);
    q_solns.push_back(q_soln);        
    return is_singular;
}

