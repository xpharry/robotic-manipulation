#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 19:44:27 2017

@author: pengxu
"""

import numpy as np
import math


# DH matrix of each coordinate
def compute_DH_matrix(d, theta, a, alpha):
    Rz_theta = np.matrix([[math.cos(theta), -math.sin(theta), 0, 0],
                          [math.sin(theta), math.cos(theta), 0, 0],
                          [0, 0, 1, 0],
                          [0, 0, 0, 1]])
    Tz_d = np.matrix([[1, 0, 0, 0],
                      [0, 1, 0, 0],
                      [0, 0, 1, d],
                      [0, 0, 0, 1]])
    Tx_a = np.matrix([[1, 0, 0, a],
                      [0, 1, 0, 0],
                      [0, 0, 1, 0],
                      [0, 0, 0, 1]])
    Rx_alpha = np.matrix([[1, 0, 0, 0],
                          [0, math.cos(alpha), -math.sin(alpha), 0],
                          [0, math.sin(alpha), math.cos(alpha), 0],
                          [0, 0, 0, 1]])
    return Rz_theta * Tz_d * Tx_a * Rx_alpha


# multiply all the DH matrice
def compute_gst_by_dh(d, theta, a, alpha):
    A = np.eye(4)
    for i in range(len(d)):
        A = A * compute_DH_matrix(d[i], theta[i], a[i], alpha[i])
    return A


def verify_se3(r):
    return r*r.transpose(), r.transpose()*r, np.linalg.det(r)


if __name__ == '__main__':
    # geometric parameters
    l0 = 13
    l1 = 14.7
    l2 = 12
    l3 = 12
    l4 = 9

    # zero configuration
    start_theta_list = np.array([0, -0.003, -0.002, 0.000, 0.000, 0.000, -1.571]).reshape(-1, 1)
    # end configuration
    end_theta_list = np.array([0, 0.122, -0.230, 1.170, -0.023, 0.750, 3.120]).reshape(-1, 1)

    # D-H Parameters
    d = [l1, 0, 0, 0, 0, 0, l4]
    theta = np.array([-math.pi/2, math.pi/2, 0, 0, 0, math.pi/2, 0]).reshape(-1, 1)
    a = [l0, 0, 0, l2, l3, 0, 0]
    alpha = [0, math.pi/2, math.pi/2, 0, 0, math.pi/2, math.pi/2]
    gst_dh = compute_gst_by_dh(d, theta+end_theta_list-start_theta_list, a, alpha)

    # turn pi/2 radius to be identical with the Tool Frame in the Exponential method
    Rz_halfpi = np.matrix([[math.cos(math.pi/2), -math.sin(math.pi/2), 0, 0],
                          [math.sin(math.pi/2), math.cos(math.pi/2), 0, 0],
                          [0, 0, 1, 0],
                          [0, 0, 0, 1]])

    gst_dh = gst_dh*Rz_halfpi
    print("gst by dh = \n", gst_dh)

    r = gst_exp[0:3, 0:3]
    r1, r2, r3 = verify_se3(r)
    print("R'R = ", r1)
    print("RR' = ", r2)
    print("|R| = ", r3)

