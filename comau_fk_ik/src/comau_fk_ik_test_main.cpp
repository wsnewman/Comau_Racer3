// comau_fk_ik_test_main.cpp
// wsn, Feb 2017
// test function for comau kinematics library

#include <comau_fk_ik/comau_fk_ik.h>

int main(int argc, char **argv) {
    ros::init(argc, argv, "comau_kinematics_test_main");
    Vectorq6x1 q_in,q_comau,q_DH;
    //q_comau << 0,-M_PI/2.0,0,0,0,0;
    for (int i=0;i<6;i++) q_in[i] = 0.0; //q_sign_corrections[i]*q_comau[i]+DH_q_offsets[i];
    q_in[1] = -M_PI/2.0; //point humerus "up", as in Comau calibration pose
    //q_in[0] = 1.0;
    
    Comau_fwd_solver comau_fwd_solver;
    Comau_IK_solver comau_IK_solver;
    
    //ROS_INFO("fwd kin at comau angles all = 0: ");
    std::cout << "q_DH: " << q_in.transpose() << std::endl;  
    q_comau = comau_IK_solver.q_comau_of_q_DH(q_in);
    
    Eigen::Affine3d A_comau_q0 = comau_fwd_solver.fwd_kin_solve_of_q_DH(q_in); //fwd_kin_solve

    Eigen::Matrix3d R_flange = A_comau_q0.linear();
    Eigen::Vector3d p_flange = A_comau_q0.translation();
    std::cout<< "p_flange = "<<p_flange.transpose()<<std::endl;
    std::cout<< "R_Flange: "<<std::endl;
    std::cout<<R_flange<<std::endl;
    
    Eigen::Matrix4d A_wrist;
    A_wrist = comau_fwd_solver.get_wrist_frame();
    //Eigen::Vector3d p_wrist = A_wrist.translation();
    std::cout<<"A_wrist: "<<std::endl;
    std::cout<<A_wrist<<std::endl;

    //return 0; // DEBUG
    ROS_INFO("==== Test for comau kinematics solver ===="); // << std::endl;
    int ans = 1;
    bool reachable_proposition;
    while (ans) {
        //pick legal, random values of q:
        double rval;
        for (int i = 0; i < 6; i++) {
            rval = ((double) rand()) / RAND_MAX;
            //ROS_INFO("randval: %f",rval);
            // assign random, legal values to joint coords;
            
            //DEBUG!!
            q_comau[i] = q_lower_limits[i] + (q_upper_limits[i]-q_lower_limits[i])*rval;
        }
        reachable_proposition = comau_fwd_solver.fit_joints_to_range(q_comau); //should not be needed
        if (reachable_proposition) {
            q_DH= comau_fwd_solver.q_DH_of_q_comau(q_comau);
            Eigen::Affine3d A_fwd_DH;
            A_fwd_DH = comau_fwd_solver.fwd_kin_solve_of_q_DH(q_DH); //fwd_kin_solve
            std::cout << "q_DH for fwd kin: " << q_DH.transpose() << std::endl;
            std::cout<<"A_fwd_DH: "<<std::endl;
            //DEBUG xxx
            //std::cout<<A_fwd_DH<<std::endl;
            Eigen::Matrix3d R_flange = A_fwd_DH.linear();
            Eigen::Matrix4d A_wrist;
            Eigen::Matrix3d R_hand;
            A_wrist = comau_fwd_solver.get_wrist_frame();

            int nsolns = comau_IK_solver.ik_solve_DH_angs(A_fwd_DH);
            std::cout << "number of IK solutions: " << nsolns << std::endl;

            std::vector<Vectorq6x1> q6dof_solns;
            comau_IK_solver.get_DH_solns(q6dof_solns); //q_DH space
            nsolns = q6dof_solns.size();
            double q_err;
            int i_min = -1;
            std::cout << "found " << nsolns << " solutions:" << std::endl;
            for (int i = 0; i < nsolns; i++) {
                Vectorq6x1 q_soln = q6dof_solns[i];
                comau_fwd_solver.fit_joints_to_range(q_soln);
                std::cout << q_soln.transpose() << std::endl;
                q6dof_solns[i] = q_soln;
                q_err =(q_DH-q_soln).norm();//fabs(q_in[0] - q_soln[0]) + fabs(q_in[1] - q_soln[1]) + fabs(q_in[2] - q_soln[2]);
                if (q_err < 0.000001) {
                    //std::cout<<"precise fit for soln "<<i<<std::endl;
                    i_min = i;
                }

                //std::cout<< "q_err: "<<q_err<<std::endl;
            }
            std::cout << "precise fit for soln " << i_min << std::endl <<std::endl;
            std::cout << "des fwd kin wrist point: " << A_wrist(0, 3) << ", " << A_wrist(1, 3) << ", " << A_wrist(2, 3) << std::endl;
            std::cout << "fwd kin wrist points from these solutions:" << std::endl;
            for (int i = 0; i < nsolns; i++) {
                A_fwd_DH = comau_fwd_solver.fwd_kin_solve_of_q_DH(q6dof_solns[i]);
                A_wrist = comau_fwd_solver.get_wrist_frame();
                std::cout << "fwd kin wrist point: " << A_wrist(0, 3) << ", " << A_wrist(1, 3) << ", " << A_wrist(2, 3) << std::endl;
            }
            ROS_INFO("fwd kin flange origins: ");
            for (int i = 0; i < nsolns; i++) {
                A_fwd_DH = comau_fwd_solver.fwd_kin_solve_of_q_DH(q6dof_solns[i]);
                std::cout << "fwd kin wrist point: " << A_fwd_DH(0, 3) << ", " << A_fwd_DH(1, 3) << ", " << A_fwd_DH(2, 3) << std::endl;
            }     
            Eigen::Vector3d bvec,nvec,tvec;
            Eigen::Matrix3d R;
            ROS_INFO("fwd kin flange b_vecs: ");
            for (int i = 0; i < nsolns; i++) {
                A_fwd_DH = comau_fwd_solver.fwd_kin_solve_of_q_DH(q6dof_solns[i]);
                R = A_fwd_DH.linear();
                bvec = R.col(2);
                std::cout << bvec.transpose()<< std::endl;
            }   
            ROS_INFO("fwd kin flange n_vecs");
            for (int i = 0; i < nsolns; i++) {
                A_fwd_DH = comau_fwd_solver.fwd_kin_solve_of_q_DH(q6dof_solns[i]);
                R = A_fwd_DH.linear();
                nvec = R.col(0);
                std::cout << nvec.transpose()<< std::endl;
            }             
            std::cout << "enter 1 to continue, 0 to stop: ";
            std::cin >> ans;
        }

    }
    return 0;
}
