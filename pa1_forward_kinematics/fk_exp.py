#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 19:44:27 2017

@author: pengxu
"""

import numpy as np
import math


# compute the (x)^, namely the skew symmetric matrix
def skew(x):
    return np.matrix([[0, -x[2], x[1]],
                      [x[2], 0, -x[0]],
                      [-x[1], x[0], 0]])


# compute the rotation matrix by w and theta
def compute_exp_w(w, theta):
    w_hat = skew(w)
    return np.eye(3) + w_hat * math.sin(theta) + w_hat * w_hat * (1 - math.cos(theta))


# compute transformation matrix
def compute_exp_epsilon(exp_w, q):
    p = np.matmul(np.eye(3) - exp_w, q)
    return np.vstack((np.hstack((exp_w, p)), np.array([0, 0, 0, 1])))


# compute the final gst expression except for gst0
def compute_gst_by_exp(w, q, theta):
    gst = np.eye(4)
    for i in range(len(w)):
        exp_ep = compute_exp_w(w[i], theta[i])
        gst = np.matmul(gst, compute_exp_epsilon(exp_ep, q[i]))
    return gst


def verify_se3(r):
    return r*r.transpose(), r.transpose()*r, np.linalg.det(r)


if __name__ == '__main__':
    # geometric parameters
    l0 = 13
    l1 = 14.7
    l2 = 12
    l3 = 12
    l4 = 9

    gst_0 = np.matrix([[1, 0, 0, l2+l3+l4],
                       [0, 1, 0, -l0],
                       [0, 0, 1, l1],
                       [0, 0, 0, 1]])

    w0 = np.array([0, 0, 1]).reshape(3, 1)
    w1 = np.array([0, -1, 0]).reshape(3, 1)
    w2 = np.array([0, 0, -1]).reshape(3, 1)
    w3 = np.array([0, 0, -1]).reshape(3, 1)
    w4 = np.array([0, 0, -1]).reshape(3, 1)
    w5 = np.array([1, 0, 0]).reshape(3, 1)

    q0 = np.array([0, -l0, l1]).reshape(3, 1)
    q1 = np.array([0, -l0, l1]).reshape(3, 1)
    q2 = np.array([0, -l0, l1]).reshape(3, 1)
    q3 = np.array([l2, -l0, l1]).reshape(3, 1)
    q4 = np.array([l2+l3, -l0, l1]).reshape(3, 1)
    q5 = np.array([l2+l3+l4, -l0, l1]).reshape(3, 1)

    w_list = [w0, w1, w2, w3, w4, w5]
    q_list = [q0, q1, q2, q3, q4, q5]
    start_theta_list = np.array([-0.003, -0.002, 0.000, 0.000, 0.000, -1.571]).reshape(-1, 1)
    end_theta_list = np.array([0.122, -0.230, 1.170, -0.023, 0.750, 3.120]).reshape(-1, 1)

    gst_exp = compute_gst_by_exp(w_list, q_list, theta=end_theta_list-start_theta_list) * gst_0
    print("gst by exp = \n", gst_exp)

    r = gst_exp[0:3, 0:3]
    r1, r2, r3 = verify_se3(r)
    print("R'R = ", r1)
    print("RR' = ", r2)
    print("|R| = ", r3)



