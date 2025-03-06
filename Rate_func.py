import numpy as np
from scipy.special import comb  # For binomial coefficient

def Rate_Factory(q_BSM, N, delta_t, L_0_in):
    L_o = (np.sqrt(3) / 3) * L_0_in

    etha_c = 0.95
    L_att = 20  # in km
    q_link_values = 0.5 * etha_c**2 * np.exp(-L_o / L_att)
    
    ExpT = np.zeros_like(q_link_values)  # To store E[T] for each q_link

    # Compute E[T] for each q_link
    for idx, q_link in enumerate(q_link_values):
        sum_term = 0

        # Compute the summation term
        for j in range(1, N + 1):
            binomial_coeff = comb(N, j, exact=True)  # Compute binomial coefficient (N choose j)
            fraction_term = 1 / (1 - (1 - q_link)**j)
            sum_term += (-1)**(j + 1) * binomial_coeff * fraction_term

        # Compute E[T] using the given expression
        ExpT[idx] = (delta_t / q_BSM**N) * sum_term

    Rate_Overall = 1 / ExpT
    return Rate_Overall


def Rate_2D(q_BSM, N, delta_t, L_0_in, m):
    L_o = (np.sqrt(3) / 3 / 2**m) * L_0_in

    etha_c = 0.95
    L_att = 20  # in km
    q_link_values = 0.5 * etha_c**2 * np.exp(-L_o / L_att)
    ETmax = np.zeros_like(q_link_values)  # to store E[T_max] for each q_link

    # Summation cutoff and tolerance for convergence in k-sum
    Kmax = 10000
    tol = 1e-8  # if term becomes very small, we stop summing

    # Loop over q_link values
    for idx, q_link in enumerate(q_link_values):
        sum_ETmax = 0

        # Sum over time steps, k (each k corresponds to time k * delta_t)
        for k in range(1, Kmax + 1):
            # Compute the CDF of T for one system at time t = k * delta_t
            F_T_k = 0
            for u in range(1, k + 1):
                p_n1 = (1 - q_BSM**N)**(u - 1) * q_BSM**N
                exponent = k // u
                # The CDF for a single geometric (with parameter q_link)
                # is 1 - (1 - q_link)^v, so for the maximum of N independent ones:
                F_n2_max = (1 - (1 - q_link)**exponent)**N
                F_T_k += p_n1 * F_n2_max
            
            # Probability that a single system hasn't finished by time k * delta_t
            term = 1 - F_T_k**(N**m)
            sum_ETmax += term

            # Check for convergence
            if term < tol:
                break

        ETmax[idx] = delta_t * sum_ETmax

    Rate_link = q_BSM**(N * (N - 1) / 2) / ETmax
    return Rate_link
