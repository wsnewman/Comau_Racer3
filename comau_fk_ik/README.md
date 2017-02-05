# comau_fk_ik
This package contains a library for computing forward and inverse kinematics for a Comau Racer-3 robot.
Key functions are:

Forward kinematics:
Eigen::Affine3d fwd_kin_solve_q_comau(const Vectorq6x1& q_comau); //use q_vec per Comau defs for this fnc

In the above, provide a vector of joint angles (consistent w/ Comau joint-angle definitions);
Return value is an affine for the flange pose.

and inverse kinematics:
    int ik_solve_comau_angs(Eigen::Affine3d const& desired_hand_pose,std::vector<Vectorq6x1> &q_comau_solns); 

Provide desired hand pose, and get back a vector full of q_vec solutions, expressed in terms
of angles consistent with Comau definitions


## Example usage
See comau_fk_ik_test_main.cpp, which is a test function to evaluate fk/ik

## Running tests/demos
with a roscore running, run:
`rosrun  comau_fk_ik comau_fk_ik_test_main`
    
