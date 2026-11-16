"""Multicomponent NRTL: tau_ij=b_ij/RT, G_ij=exp(-alpha_ij tau_ij)
ln g_i = sum_j x_j tau_ji G_ji / sum_k x_k G_ki
         + sum_j [x_j G_ij / sum_k x_k G_kj] (tau_ij - sum_m x_m tau_mj G_mj / sum_k x_k G_kj)"""
import numpy as np
from .params import R, pair

def tau_G(comps, T):
    n = len(comps); tau = np.zeros((n, n)); G = np.ones((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                p = pair(comps[i], comps[j])
                tau[i, j] = p[2] / (R * T); G[i, j] = np.exp(-p[4] * tau[i, j])
    return tau, G

def gamma_nrtl(x, T, comps):
    x = np.asarray(x, float); tau, G = tau_G(comps, T); n = len(x)
    S = G.T @ x                 # S_j = sum_k x_k G_kj
    C = (tau * G).T @ x         # C_j = sum_m x_m tau_mj G_mj
    lng = np.zeros(n)
    for i in range(n):
        t1 = C[i] / S[i]
        t2 = sum(x[j] * G[i, j] / S[j] * (tau[i, j] - C[j] / S[j]) for j in range(n))
        lng[i] = t1 + t2
    return np.exp(lng)
