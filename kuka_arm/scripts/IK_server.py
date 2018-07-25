#!/usr/bin/env python

# Copyright (C) 2017 Udacity Inc.
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
# All Rights Reserved.
#
# Author: Grace Livingston

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from sympy import *
from math import pi

def create_mats():
    """Create the DH Params and rotation matrices"""

    # Create symbols for DH param

    # Joint angles (Theta)
    q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')

    # Link offsets
    d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')

    # Link lengths
    a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')

    # Joint twists (alpha)
    p0, p1, p2, p3, p4, p5, p6 = symbols('p0:7')

    # DH Parameters
    dh = {p0:      0, a0:      0, d1:  0.75, q1:        q1,
          p1: -pi/2., a1:   0.35, d2:     0, q2: -pi/2.+q2,
          p2:      0, a2:   1.25, d3:     0, q3:        q3,
          p3: -pi/2., a3: -0.054, d4:   1.5, q4:        q4,
          p4:  pi/2., a4:      0, d5:     0, q5:        q5,
          p5: -pi/2., a5:      0, d6:     0, q6:        q6,
          p6:      0, a6:      0, d7: 0.303, q7:         0}

    # Define Modified DH Transformation matrix
    # Function to return homogeneous transform matrix
    def TF_Mat(p, a, d, q):
        return Matrix([
            [        cos(q),       -sin(q),       0,         a],
            [ sin(q)*cos(p), cos(q)*cos(p), -sin(p), -sin(p)*d],
            [ sin(q)*sin(p), cos(q)*sin(p),  cos(p),  cos(p)*d],
            [             0,             0,       0,         1]
        ])

    # Transform matrices
    T0_1 = TF_Mat(p0, a0, d1, q1).subs(dh)
    T1_2 = TF_Mat(p1, a1, d2, q2).subs(dh)
    T2_3 = TF_Mat(p2, a2, d3, q3).subs(dh)
    T3_4 = TF_Mat(p3, a3, d4, q4).subs(dh)
    T4_5 = TF_Mat(p4, a4, d5, q5).subs(dh)
    T5_6 = TF_Mat(p5, a5, d6, q6).subs(dh)
    T6_7 = TF_Mat(p6, a6, d7, q7).subs(dh)

    # Homogenous Transform from link_0 to end link_7
    T0_7 = (T0_1 * T1_2 * T2_3 * T3_4 * T4_5 * T5_6 * T6_7)

    # Correction from URDF to DH
    R_corr = Matrix([
        [  0.0,  0.0,  1.0,  0.0],
        [  0.0, -1.0,  0.0,  0.0],
        [  1.0,  0.0,  0.0,  0.0],
        [  0.0,  0.0,  0.0,  1.0]
    ])

    # HT from 0 - 7 with URDF - DH correction.
    T0_7_corr = (T0_7 * R_corr)

    # End Effector rotation matrix
    r,p,y = symbols('r p y')

    ROT_x = Matrix([[       1,       0,       0],
                    [       0,  cos(r), -sin(r)],
                    [       0,  sin(r),  cos(r)]])
    ROT_y = Matrix([[  cos(p),       0,  sin(p)],
                    [       0,       1,       0],
                    [ -sin(p),       0,  cos(p)]])
    ROT_z = Matrix([[  cos(y), -sin(y),       0],
                    [  sin(y),  cos(y),       0],
                    [       0,       0,       1]])

    # Correction for Gazebo - DH Parameters
    ROT_corr = ROT_z.subs(y, pi) * ROT_y.subs(p, -pi / 2.)

    ROT_EE = ROT_z * ROT_y * ROT_x * ROT_corr

    return ROT_EE, T0_1, T1_2, T2_3

# IK handler service
def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:

        ROT_EE = create_mats()
        # Initialize service response
        joint_trajectory_list = []
        for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

            # px,py,pz = end-effector position
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            EE = Matrix([
                [px],
                [py],
                [pz]
            ])

            # roll, pitch, yaw = end-effector orientation
            (roll,pitch,yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x,
                 req.poses[x].orientation.y,
                 req.poses[x].orientation.z,
                 req.poses[x].orientation.w])

            ROT_EE = ROT_EE.subs({'r': roll, 'p': pitch, 'y': yaw})

            # Wrist center
            WC = EE - (0.303) * ROT_EE[:,2]

            # Compute joint angles
            # Calculate Theta1
            theta1 = atan2(WC[1],WC[0])

            # find the 3rd side of the triangle
            B = sqrt(
                pow((sqrt(WC[0]*WC[0] + WC[1]*WC[1]) - 0.35), 2) +
                pow((WC[2] - 0.75), 2)
            )

            # Triangle angles
            a = acos((-0.6875 + B*B) / (2.5*B))
            b = acos(( 3.8125 - B*B) / (3.75))
            c = acos(( 0.6875 + B*B) / (3.0*B))

            # Find Theta2 and Theta3
            theta2 = (
                pi / 2 -
                a -
                atan2(WC[2]-0.75,
                      sqrt(WC[0] * WC[0] +
                           WC[1] * WC[1]) -
                      0.35)
            )

            # 0.036 compensates for link 4 sag of -0.054m
            theta3 = pi / 2 - (b+0.036)

            # Extract rotation matrices from the transformation matrices
            # Extract rotation matrix R0_3 from transformation matrix T0_3 then substitute angles q1-3
            R0_3 = T0_1[0:3,0:3] * T1_2[0:3,0:3] * T2_3[0:3,0:3]
            R0_3 = R0_3.evalf(subs={q1: theta1, q2: theta2, q3:theta3})

            # Get R3_6 rotation matrix
            R3_6 = R0_3.transpose() * ROT_EE

            # Euler angles from rotation matrix
            theta5 = atan2(
                sqrt(R3_6[0,2] * R3_6[0,2] +
                     R3_6[2,2] * R3_6[2,2]),
                R3_6[1,2]
            )

            # use Theta5 to find best solution
            if (theta5 > pi) :
                theta4 = atan2(-R3_6[2,2],  R3_6[0,2])
                theta6 = atan2( R3_6[1,1], -R3_6[1,0])
            else:
                theta4 = atan2( R3_6[2,2], -R3_6[0,2])
                theta6 = atan2(-R3_6[1,1],  R3_6[1,0])

            # Populate response for the IK request
            #
            # In the next line replace theta1,theta2...,theta6 by your
            # joint angle variables
            joint_trajectory_point.positions = [
                theta1,
                theta2,
                theta3,
                theta4,
                theta5,
                theta6
            ]
            joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" %
                      len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)

def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
