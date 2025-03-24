import numpy as np
from scipy.special import comb  # For binomial coefficient

def Rate_Factory(q_BSM, N, delta_t, L_0_in):
    L_o = L_0_in / (2 * np.sin(np.pi / N))

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



def Rate_Decent(q_BSM, q_Fuse, N, delta_t, L_0_in, k_max):
    """
    Rate_Decent calculates the average GHZ entanglement distribution rate
    for an N-qubit system using a decentralized switching approach.

    Parameters:
      q_BSM    - Success probability of a Bell-state measurement (BSM)
      q_Fuse   - Success probability of a fusion operation
      N        - Number of qubits/nodes involved in GHZ entanglement
      delta_t  - Time duration of a single attempt (time step in s)
      L_0_in   - The final distance between neighboring nodes (in km)
      k_max    - Maximum number of time slots to consider

    Returns:
      Rate_out - Average entanglement distribution rate (inverse of E[T_max])
    """
    
    L_o = L_0_in / 2   # element-wise division if L_0_in is an array
    etha_c = 0.95      # Coupling efficiency
    L_att  = 20        # Attenuation length (in km)

    # Compute overall link success probability
    q_link_values = 0.5 * etha_c**2 * np.exp(-L_o / L_att)
    
    # Ensure q_link_values is array-like
    q_link_values = np.atleast_1d(q_link_values)

    # Initialize an array to store the expected time E[T] for each q_link
    E_Tmax = np.zeros(q_link_values.shape)

    # Loop over each q_link value
    for idx, q_link in np.ndenumerate(q_link_values):
        
        # Initialize arrays for the CDF values (using 0-indexing)
        F_n_i = np.zeros(k_max)  # CDF of a single n_i (individual system completion time)
        F_n   = np.zeros(k_max)  # CDF of the maximum completion time: n = max{n_1, ..., n_N}
        
        # Each n_i = n_BSM * n_dist, where:
        # - n_BSM follows a geometric distribution with success probability q_BSM
        # - n_dist has CDF: [1 - (1 - q_link)^v]^2 evaluated at v = floor(k/u)
        for k in range(1, k_max + 1):
            sum_u = 0.0
            for u in range(1, k + 1):
                # Probability that n_BSM = u (geometric distribution)
                p_nBSM = (1 - q_BSM)**(u - 1) * q_BSM
    
                # For fixed n_BSM = u, compute L = floor(k/u)
                L_inside = k // u  # integer division
                # CDF of n_dist at L
                p_nDist_CDF = (1 - (1 - q_link)**L_inside)**2
    
                # Add contribution to the overall CDF of n_i
                sum_u += p_nBSM * p_nDist_CDF
            F_n_i[k - 1] = sum_u            # CDF of n_i up to time slot k
            F_n[k - 1] = F_n_i[k - 1]**N     # CDF of max{n_i} using independence
        
        E_n = 0.0
        for k in range(1, k_max + 1):
            E_n += 1 - F_n[k - 1]
        
        # Multiply by the time step delta_t and account for fusion operations.
        # Each fusion step must succeed independently for all N fusions.
        E_Tmax[idx] = (E_n * delta_t) / (q_Fuse**N)
    
    Rate_out = 1 / E_Tmax  # Inverse of average time gives the entanglement rate
    
    # If the output is a single value array, return a scalar
    if Rate_out.size == 1:
        return Rate_out.item()
    return Rate_out
