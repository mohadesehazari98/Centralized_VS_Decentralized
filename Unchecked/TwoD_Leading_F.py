"""
Ref: This code has been motivated by (Leading_F.py),
the original authorship is preserved to respect the initial authors.
Mohadeseh Azari / Department of Informatics and Networked Systems / 
School of Computing and Information / University of Pittsburgh / 
Pittsburgh,PA / moa125@pitt.edu'
"""
from Leading_F import fidelity_to_depolarizing_prob, _check_set
import itertools


def _g_function(num_remote_nodes, GHZ_success_prob, subset, prob):
    """This _g_function is the same as Leading-order expression for the function G in Appendix D of the paper
    "Analysis of Multipartite Entanglement Distribution using a Central Quantum-Network Node"
    Guus Avis, Filip Rozpędek and Stephanie Wehner.

    Parameters
    ----------
    num_remote_nodes : int
        Number of remote nodes that are not yet in a maximally mixed state.
    GHZ_success_prob : float
        Probability per attempt that entanglement distribution across Parent's node succeeds. 
        IOW success probability of generating a Parent == Rate .
    subset : set of int
        A subset of {1, 2, ..., num_remote_nodes} = W.
    prob : float
        Number between 0 and 1, typically a depolarizing probability (1 - p_mem ** N).
    """
    _check_set(subset, num_remote_nodes)
    if num_remote_nodes == 0:
        return 0
    value = 1
    for i in range(1, num_remote_nodes + 1):
        value *= (num_remote_nodes + 1 - i) * GHZ_success_prob
        value /= sum([1 for j in subset if j < i]) * prob + (num_remote_nodes + 1 - i) * GHZ_success_prob
    return value


def G_function(num_remote_nodes, set_all_integers, GHZ_success_prob, prob, link_bsm_depolar):
    """This G_function is defined in the notes.

    Parameters
    ----------
    num_remote_nodes : int
        Number of remote nodes that are not yet in a maximally mixed state.
    set_all_intergers:
        As the total number of end nodes that wish to share a GHZ state among them 
    GHZ_success_prob : float
        Probability per attempt that entanglement distribution across Parent's node succeeds. 
        IOW success probability of generating a Parent == Rate .
    prob : float
        Number between 0 and 1, typically a depolarizing probability (1 - p_mem ** N).
    link_bsm_depolar : float 
        The probability of all (N-1) entanglement swap success in every Parent 
    """
    g_func = _g_function
    fidelity = 0
    for i in range(num_remote_nodes+1):
        i_even = i % 2 == 0
        if not i_even and i != num_remote_nodes:
            continue
        prefactor = (1 / 2) ** num_remote_nodes if i_even else 0
        if i == num_remote_nodes:
            prefactor += 1 / 2
        subsets_length_i = itertools.combinations(set_all_integers, i)
        gs = 0
        for subset in subsets_length_i:
            g = g_func(num_remote_nodes=num_remote_nodes,
                       GHZ_success_prob=GHZ_success_prob,
                       subset=set(subset),
                       prob=prob)
            gs += g
        fidelity += prefactor * link_bsm_depolar ** i * gs
    return fidelity


def TwoD_network_fidelity(num_remote_nodes, GHZ_success_prob, mem_depolar_prob, bsm_depolar_prob, ghz_fidelity):
    """Analytical results for the fidelity achieved with the 2D repeater protocol (leading-order expression).

    The leading-order expression is to leading order in the link success probability (`link_success_prob`)
    and the probability that a qubit undergoes a depolarizing error in memory while performing a single attempt
    at Bell-state distribution (`mem_depolar_prob`).

    Parameters
    ----------
    num_remote_nodes : int
        Number of remote nodes that attempt to share a GHZ state through the GHZ-factory protocol.
    GHZ_success_prob : float
        Probability per attempt that entanglement distribution across Parent's node succeeds. 
        IOW success probability of generating a Parent == Rate .
    ghz_fidelity : float
        Fidelity of the GHZ state that each Parent possess.
    mem_depolar_prob : float
        Depolarizing noise in quantum memory per time step for both the GHZ factory and end nodes.
    bsm_depolar_prob : float
        Depolarizing probability on each qubit participating in a Bell-state measurement.

    """

    num_remote_nodes = int(num_remote_nodes)  # in case it's passed as a float
    ghz_depolarizing_prob = fidelity_to_depolarizing_prob(num_qubits=num_remote_nodes, fidelity=ghz_fidelity)
    link_bsm_depolar = (1 - bsm_depolar_prob) ** (num_remote_nodes - 1)
    prob = 1 - (1 - mem_depolar_prob) ** (num_remote_nodes)
    set_all_integers = set(range(1, num_remote_nodes + 1))

    fidelity = 0
    for i in range(num_remote_nodes + 1):
        subsets_length_i = itertools.combinations(set_all_integers, i)
        gs = 0
        if i == 0:
            gs = 1
        else:
            for subset in subsets_length_i:
                g = G_function(num_remote_nodes=int(i), 
                            set_all_integers=set(range(1, int(i) + 1)), 
                            GHZ_success_prob=GHZ_success_prob, 
                            prob=prob, 
                            link_bsm_depolar=link_bsm_depolar)
                gs += g
        fidelity += ((1/2 - ghz_depolarizing_prob/2) ** (num_remote_nodes - i)) * (ghz_depolarizing_prob ** i) * gs 
    return fidelity

