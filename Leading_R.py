import numpy as np
import matplotlib.pyplot as plt
import math

def factory_waiting_time(num_remote_nodes, link_success_prob, bsm_success_prob):
    """Leading-order expression of the average waiting time to deliver a single GHZ state with the GHZ-factory protocol.

    The expression is to leading order in the link success probability.
    It thus becomes exact in the limit when the link success probability goes to zero.

    Parameters
    ----------
    num_remote_nodes : int
        Number of remote nodes that attempt to share a GHZ state through the GHZ-factory protocol.
    link_success_prob : float
        Probability per attempt that entanglement distribution between a remote node and the factory succeeds.
    bsm_success_prob : float
        Success probability of Bell-state measurements.

    """

    num_remote_nodes = int(num_remote_nodes)  # in case it's passed as a float
    # find harmonic number H_N where N = num_remote_nodes
    

    distribution_time = n_all(num_remote_nodes, link_success_prob) / (bsm_success_prob ** num_remote_nodes)

    return 1 / distribution_time


def n_all(num_remote_nodes, link_success_prob):
    total_sum = 0
    for j in range(1, num_remote_nodes + 1):
        term = ((-1) ** (j + 1)) * math.comb(num_remote_nodes, j) * (1 / (1 - (1 - link_success_prob) ** j))
        total_sum += term
    return total_sum


num_remote_nodes=5
link_success_prob = np.linspace(1e-10, 1, 100000)  # Start from a small number close to 0
bsm_success_prob = 1

rates = [factory_waiting_time(num_remote_nodes, q_link, bsm_success_prob) for q_link in link_success_prob]

# Plotting
plt.figure(figsize=(8, 5))
plt.plot(link_success_prob, rates, label="Rate vs Link Success Prob", color='b')
plt.xlabel("Link Success Probability")
plt.ylabel("Rate (R)")
plt.ylim((0, 1))
plt.title("Rate as a Function of Link Success Probability")
plt.legend()
plt.grid(True)
plt.show()