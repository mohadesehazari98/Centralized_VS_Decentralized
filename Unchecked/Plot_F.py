from Leading_F import factory_fidelity
from TwoD_Leading_F import TwoD_network_fidelity
from Rate_func import Rate_Factory
import numpy as np
import matplotlib.pyplot as plt

# Define parameters
L_0_in = np.linspace(1e-3, 500, 100000)  # Start from a small number close to 0
mem_depolar_prob = 1e-2
link_depolar_prob = 0
bsm_depolar_prob = 0
bsm_success_prob = 1 - 1e-4
ghz_fidelity = 1

# Define different colors for plotting
colors = ['b', 'r']

# Plotting
plt.figure(figsize=(8, 5))

for idx, num_remote_nodes in enumerate([3, 4]):
    # Compute GHZ success probability
    GHZ_success_prob = Rate_Factory(q_BSM=bsm_success_prob, N=num_remote_nodes, 
                                    delta_t=1, L_0_in=L_0_in)

    # Compute F_target fidelity
    F_target = factory_fidelity(num_remote_nodes, L_0_in, mem_depolar_prob, link_depolar_prob,
                                bsm_depolar_prob, ghz_fidelity)

    # Compute F_TwoDim fidelity
    F_TwoDim = TwoD_network_fidelity(num_remote_nodes, GHZ_success_prob, mem_depolar_prob, 
                                     bsm_depolar_prob, F_target)

    # Plot F_target and F_TwoDim for each num_remote_nodes
    plt.plot(L_0_in, F_target, label=f"F_target (N={num_remote_nodes})", color=colors[idx])
    plt.plot(L_0_in, F_TwoDim, label=f"F_TwoDim (N={num_remote_nodes})", color=colors[idx], linestyle='dashed')

plt.xlabel("Neighboring Distance (km)")
plt.ylabel("Fidelity")
plt.ylim((0, 1))
plt.title("Fidelity of the child and its parent for 3-GHZ and 4-GHZ")
plt.legend()
plt.grid(True)
plt.show()
