#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thur Oct 12 13:54:27 2017

@author: Peng Xu
"""

import numpy as np
from numpy.linalg import inv, norm
import math

# const
l0 = 13
l1 = 14.7
l2 = 12
l3 = 12
l4 = 9

eps = 1e-2


# compute the (x)^, namely the skew symmetric matrix
def skew(x):
    return np.matrix([[0, -x[2], x[1]],
                      [x[2], 0, -x[0]],
                      [-x[1], x[0], 0]])


# compute the rotation matrix by w and theta
def compute_exp_w(w, theta):
    w_hat = skew(w)
    return np.matrix(np.eye(3)) + w_hat * math.sin(theta) + w_hat * w_hat * (1 - math.cos(theta))


# compute transformation matrix
def compute_exp_epsilon(w, theta, q):
    exp_w = compute_exp_w(w, theta)
    p = np.matmul(np.eye(3) - exp_w, q)
    return np.matrix(np.vstack((np.hstack((exp_w, p)), np.array([0, 0, 0, 1]))))


def subproblem_1(w, r, p, q):
    r = np.matrix(r.tolist()[:3]).reshape(-1, 1)
    p = np.matrix(p.tolist()[:3]).reshape(-1, 1)
    q = np.matrix(q.tolist()[:3]).reshape(-1, 1)

    u = p - r
    v = q - r
    u_prime = u - w * w.transpose() * u
    v_prime = v - w * w.transpose() * v

    if abs(w.transpose()*u - w.transpose()*v) < eps \
            and abs(u_prime.transpose()*u_prime - v_prime.transpose()*v_prime) < eps:
        return np.asscalar(np.arctan2(w.transpose()*np.cross(u_prime, v_prime, axis=0),
                                      u_prime.transpose()*v_prime))
    else:
        print("********** Debug **********")
        print(abs(w.transpose()*u - w.transpose()*v))
        print(abs(u_prime.transpose()*u_prime - v_prime.transpose()*v_prime))
        return np.asscalar(np.arctan2(w.transpose()*np.cross(u_prime, v_prime, axis=0),
                                      u_prime.transpose()*v_prime))


def subproblem_2(w1, w2, r, p, q):
    r = np.matrix(r.tolist()[:3]).reshape(-1, 1)
    p = np.matrix(p.tolist()[:3]).reshape(-1, 1)
    q = np.matrix(q.tolist()[:3]).reshape(-1, 1)

    u = p - r
    v = q - r
    temp = np.square(w1.transpose()*w2) - 1
    if abs(temp) < eps:
        return [], []
    alpha = ((w1.transpose()*w2)*w2.transpose()*u - w1.transpose()*v) / temp
    beta = ((w1.transpose()*w2)*w1.transpose()*v - w2.transpose()*u) / temp
    gama_square = (u.transpose()*u - np.square(alpha) - np.square(beta) - 2*alpha*beta*w1.transpose()*w2) \
                  / np.matmul(np.cross(w1, w2, axis=0).transpose(), np.cross(w1, w2, axis=0))
    alpha = np.asscalar(alpha)
    beta = np.asscalar(beta)
    gama_square = np.asscalar(gama_square)

    if abs(gama_square) < eps:
        z = alpha*w1 + beta*w2
        c = z + r
        theta1 = subproblem_1(-w1, r, q, c)
        theta2 = subproblem_1(w2, r, p, c)
        return [theta1], [theta2]
    elif gama_square > 0:
        gama = np.sqrt(gama_square)
        # z1
        z1 = alpha*w1 + beta*w2 + gama*(np.cross(w1, w2, axis=0))
        c1 = z1 + r
        theta11 = subproblem_1(-w1, r, q, c1)
        theta21 = subproblem_1(w2, r, p, c1)
        # z2
        z2 = alpha*w1 + beta*w2 - gama*(np.cross(w1, w2, axis=0))
        c2 = z2 + r
        theta12 = subproblem_1(-w1, r, q, c2)
        theta22 = subproblem_1(w2, r, p, c2)
        return [theta11, theta12], [theta21, theta22]
    else:
        return [], []


def subproblem_3(w, r, p, q, delta):
    r = np.matrix(r.tolist()[:3]).reshape(-1, 1)
    p = np.matrix(p.tolist()[:3]).reshape(-1, 1)
    q = np.matrix(q.tolist()[:3]).reshape(-1, 1)

    u = p - r
    v = q - r
    u_prime = u - w * w.transpose() * u
    v_prime = v - w * w.transpose() * v
    delta_prime_square = delta**2 - (w.transpose()*(p-q))**2
    theta_0 = math.atan2(w.transpose() * np.cross(u_prime, v_prime, axis=0), u_prime.transpose()*v_prime)

    d_u_prime = np.sqrt(u_prime.transpose() * u_prime)
    d_v_prime = np.sqrt(v_prime.transpose() * v_prime)
    d_delta_prime = np.sqrt(delta_prime_square)

    if d_u_prime + d_v_prime < d_delta_prime or d_u_prime < eps or d_v_prime < eps:
        return []

    if abs(d_u_prime + d_v_prime - d_delta_prime) < eps:
        return [theta_0]

    offset = math.acos((u_prime.transpose()*u_prime + v_prime.transpose()*v_prime - delta_prime_square)
                       / (2*d_u_prime*d_v_prime))
    return [theta_0+offset, theta_0-offset]


# test functions
def test():
    w1 = np.matrix([1/math.sqrt(8), math.sqrt(3/8), 1/math.sqrt(2)]).reshape(-1, 1)
    w2 = np.matrix([-1/math.sqrt(3), 0, math.sqrt(2/3)]).reshape(-1, 1)
    r = np.matrix([1, 2, -3]).reshape(-1, 1)
    p = np.matrix([5, -7, 12, 1]).reshape(-1, 1)
    q1 = np.matrix([17.3037, -3.1128, 2.48175]).reshape(-1, 1)
    q2 = np.matrix([14.6711, -9.34914, -0.490328]).reshape(-1, 1)
    q3 = np.matrix([5, -3, -12]).reshape(-1, 1)
    delta = 19.0031

    print("\nTest Results:\n")
    theta_sol = subproblem_1(w1, r, p, q1)
    print("theta_sol1 = ", theta_sol)
    theta_sol = subproblem_2(w1, w2, r, p, q2)
    print("theta_sol2 = ", theta_sol)
    theta_sol = subproblem_3(w1, r, p, q3, delta)
    print("theta_sol3 = ", theta_sol)


def main():
    # g_d and g_init
    g_d = np.matrix([[-0.0256, -0.4496, -0.8929, -4.197],
                     [0.9758, 0.1829, -0.1200, 15.369],
                     [0.2173, -0.8743, 0.4340, 13.931],
                     [0, 0, 0, 1]])
    g_0 = np.matrix([[1, 0, 0, l2+l3+l4],
                     [0, 1, 0, -l0],
                     [0, 0, 1, l1],
                     [0, 0, 0, 1]])
    g1 = g_d * inv(g_0)
    g2 = g_0 * inv(g_d)

    theta_init = [-0.003, -0.002, 0.000, 0.000, 0.000, -1.571]

    # axes
    w1 = np.matrix([0, 0, 1]).reshape(-1, 1)
    w2 = np.matrix([0, -1, 0]).reshape(-1, 1)
    w3 = np.matrix([0, 0, -1]).reshape(-1, 1)
    w4 = np.matrix([0, 0, -1]).reshape(-1, 1)
    w5 = np.matrix([0, 0, -1]).reshape(-1, 1)
    w6 = np.matrix([1, 0, 0]).reshape(-1, 1)

    q1 = np.matrix([0, -l0, l1]).reshape(-1, 1)
    q2 = np.matrix([0, -l0, l1]).reshape(-1, 1)
    q3 = np.matrix([0, -l0, l1]).reshape(-1, 1)
    q4 = np.matrix([l2, -l0, l1]).reshape(-1, 1)
    q5 = np.matrix([l2+l3, -l0, l1]).reshape(-1, 1)
    q6 = np.matrix([l2+l3, -l0, l1]).reshape(-1, 1)

    # define reference points
    p1 = np.matrix([0, -l0, l1, 1]).reshape(-1, 1)  # intersection point of eta1, eta2 and eta3
    p2 = np.matrix([l2, -l0, l1, 1]).reshape(-1, 1)  # a point on eta4
    p3 = np.matrix([l2+l3, -l0, l1, 1]).reshape(-1, 1)  # intersection point of eta5 and eta6

    theta_sols = []

    # Step 1: SP3
    delta = g2*p1 - p3
    delta = norm(delta)
    theta4_sol = subproblem_3(-w4, p2, p1, p3, delta)
    for i in range(len(theta4_sol)):
        theta4 = theta4_sol[i]
        # print("theta4 = {}".format(theta4))

        # Step 2: SP2
        e_4_inv = inv(compute_exp_epsilon(w4, theta4, q4))
        p4 = e_4_inv * p1
        theta6_sol, theta5_sol = subproblem_2(-w6, -w5, p3, p4, g2*p1)
        for j in range(len(theta5_sol)):
            theta5 = theta5_sol[j]
            theta6 = theta6_sol[j]
            # print("theta5 = {}, theta6 = {}".format(theta5, theta6))

            # Step 3: SP2
            p5 = np.array([0, -l0, 0, 1]).reshape(-1, 1)
            e_5_inv = inv(compute_exp_epsilon(w5, theta5, q5))
            e_6_inv = inv(compute_exp_epsilon(w6, theta6, q6))
            g3 = g1 * e_6_inv * e_5_inv * e_4_inv
            theta1_sol, theta2_sol = subproblem_2(w1, w2, p1, p5, g3*p5)
            for k in range(len(theta1_sol)):
                theta1 = theta1_sol[k]
                theta2 = theta2_sol[k]
                # print("theta1 = {}, theta2 = {}".format(theta1, theta2))

                # Step 4:
                p6 = np.copy(p2)
                e_1_inv = inv(compute_exp_epsilon(w1, theta1, q1))
                e_2_inv = inv(compute_exp_epsilon(w2, theta2, q2))
                g4 = e_2_inv * e_1_inv * g3
                theta3 = subproblem_1(w3, p1, p6, g4*p6)
                # print("theta3 = {}".format(theta3))

                # combine results
                theta_end = [theta1, theta2, theta3, theta4, theta5, theta6]
                theta_sol = list(np.array(theta_end) - np.array(theta_init))
                theta_sols.append(theta_sol)

    print("\nFinal Results:\n")
    for i in range(len(theta_sols)):
        print("Solution ", i+1, ":")
        print(theta_sols[i])


if __name__ == '__main__':
    print("\n********** Test **********")
    test()
    print("\n********** Main **********")
    main()

