"""Multicomponent Wilson: Lambda_ij = (Vj/Vi) exp(-a_ij/RT)
ln g_i = 1 - ln(sum_j x_j L_ij) - sum_k x_k L_ki / sum_j x_j L_kj"""
import numpy as np
from .params import COMPONENTS, R, pair

def lambda_matrix(comps, T):
    n = len(comps); L = np.ones((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                a_ij = pair(comps[i], comps[j])[0]
                Vi, Vj = COMPONENTS[comps[i]]["V"], COMPONENTS[comps[j]]["V"]
                L[i, j] = Vj / Vi * np.exp(-a_ij / (R * T))
    return L

def gamma_wilson(x, T, comps):
    x = np.asarray(x, float); L = lambda_matrix(comps, T)
    s = L @ x                       # s_k = sum_j x_j L_kj
    lng = 1 - np.log(s) - (L.T @ (x / s))   # sum_k x_k L_ki / s_k
    return np.exp(lng)
