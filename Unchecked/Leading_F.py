"""
Ref: This code has been modified, and the original authorship is preserved to respect the initial authors.
----
This code can be used to evaluate the analytical results for the fidelity of GHZ-state distribution
on a symmetric star-shaped network using a factory node that are presented in the paper
"Analysis of Multipartite Entanglement Distribution using a Central Quantum-Network Node"
by Guus Avis, Filip Rozpędek and Stephanie Wehner.
It includes leading-order results for the fidelity.

We note that a different convention for the depolarizing channel is used in this file than in the paper.
In the paper, a depolarizing parameter with param p has the action on one qubit
rho -> p rho + (1 - p) Id_2 / 2.
Here, the convention that is used is
rho -> (1 - p) rho + p Id_2 / 2.
This makes it easier to work to leading order in noise parameters.

Because of the difference in convention, there is the following mapping between parameter in the paper and parameters
used here:
link_depolar_prob = 1 - p_link
bsm_depolar_prob = 1 - p_BSM
mem_depolar_prob = 1 - p_mem

Additionally, we here don't use the parameter p_GHZ for the noise in the GHZ state as it is locally prepared.
Instead, we use the fidelity of the prepared state,
ghz_fidelity = 1 - (2 ** num_qubits - 1) / 2 ** num_qubits * p_ghz
(see also documentation of the function `fidelity_to_depolarizing_prob()`).


"""
import itertools
import numpy as np


def _g_function(num_remote_nodes, link_success_prob, subset, prob):
    """Leading-order expression for the function G in Appendix D of the paper.

    For definition and derivation, see Appendix D of the paper
    "Analysis of Multipartite Entanglement Distribution using a Central Quantum-Network Node"
    Guus Avis, Filip Rozpędek and Stephanie Wehner.

    Parameters
    ----------
    num_remote_nodes : int
        Number of remote nodes that attempt to share a GHZ state through the GHZ-factory protocol.
    link_success_prob : float
        Probability per attempt that entanglement distribution between a remote node and the factory succeeds.
    subset : set of int
        A subset of {1, 2, ..., num_remote_nodes} = U.
    prob : float
        Number between 0 and 1, typically a depolarizing probability (1 - p_mem ** 2).

    """
    _check_set(subset, num_remote_nodes)
    value = 1
    for i in range(1, num_remote_nodes + 1):
        value *= (num_remote_nodes + 1 - i) * link_success_prob
        value /= sum([1 for j in subset if j < i]) * prob + (num_remote_nodes + 1 - i) * link_success_prob
    return value


def factory_fidelity(num_remote_nodes, L_0_in, mem_depolar_prob, link_depolar_prob,
                     bsm_depolar_prob, ghz_fidelity):
    """Analytical results for the fidelity achieved with the GHZ-factory protocol (leading-order expression or bound).

    Parameter `bound` can be used to toggle between a strict lower bound on the fidelity and a leading-order expression.

    The leading-order expression is to leading order in the link success probability (`link_success_prob`)
    and the probability that a qubit undergoes a depolarizing error in memory while performing a single attempt
    at Bell-state distribution (`mem_depolar_prob`).

    Parameters
    ----------
    num_remote_nodes : int
        Number of remote nodes that attempt to share a GHZ state through the GHZ-factory protocol.
    link_success_prob : float
        Probability per attempt that entanglement distribution between a remote node and the factory succeeds.
    ghz_fidelity : float
        Fidelity of the GHZ state that is created at the GHZ factory.
    mem_depolar_prob : float
        Depolarizing noise in quantum memory per time step for both the GHZ factory and end nodes.
    link_depolar_prob : float
        Depolarizing noise on EPR pairs delivered at entanglement generation.
        The state that is delivered is assumed to be (1 - link_depolar_prob)|Phi+><Phi+| + link_depolar_prob * Id / 2
    bsm_depolar_prob : float
        Depolarizing probability on each qubit participating in a Bell-state measurement.
    bound : bool (optional)
        If True, a lower bound is calculated.
        If False, a leading-order expression is calculated.

    """

    L_o = (np.sqrt(3) / 3) * L_0_in
    etha_c = 0.95
    L_att = 20  # in km
    link_success_prob = 0.5 * etha_c**2 * np.exp(-L_o / L_att)
    g_func = _g_function

    num_remote_nodes = int(num_remote_nodes)  # in case it's passed as a float
    ghz_depolarizing_prob = fidelity_to_depolarizing_prob(num_qubits=num_remote_nodes, fidelity=ghz_fidelity)
    link_bsm_depolar = (1 - link_depolar_prob) * (1 - bsm_depolar_prob) ** 2
    prob = 2 * mem_depolar_prob - mem_depolar_prob ** 2
    set_all_integers = set(range(1, num_remote_nodes + 1))

    fidelity = 0
    for i in range(num_remote_nodes + 1):
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
                       link_success_prob=link_success_prob,
                       subset=set(subset),
                       prob=prob)
            gs += g
        fidelity += prefactor * link_bsm_depolar ** i * gs

    fidelity *= 1 - ghz_depolarizing_prob

    fidelity += ghz_depolarizing_prob / 2 ** num_remote_nodes

    return fidelity


def fidelity_to_depolarizing_prob(fidelity, num_qubits):
    """Depolarizing probability needed to reach a certain fidelity by depolarizing a pure state.

    Depolarizing channel: rho -> rho' = (1 - depolarizing_prob) rho + depolarizing_prob * Id / 2 ** num_qubits.
    This function returns depolarizing_prob needed for rho' to have a specific fidelity to rho in case rho is pure.
    The fidelity of rho to rho is 1, and the fidelity of rho to Id is also 1 (because rho is pure).
    Therefore, the fidelity of rho' to rho is (1 - depolarizing_prob) + depolarizing_prob / 2 ** num_qubits.
    This can be rewritten as fidelity = 1 - (2 ** num_qubits - 1) / 2 ** num_qubits * depolarizing_prob.
    This can be inverted to give depolarizing_probability = 2 ** num_qubits / (2 ** num_qubits - 1) * (1 - fidelity).

    Parameters
    ----------
    fidelity : float
        Fidelity (squared definition) of the depolarized state to the state before the depolarizing channel.
    num_qubits : float
        Number of qubits in the state that undergoes a depolarizing channel.

    Returns
    -------
    depolarizing_prob : float
        Fidelity of the depolarized state to the state before depolarization.

    """
    return 2 ** num_qubits / (2 ** num_qubits - 1) * (1 - fidelity)


def _check_set(subset, num_remote_nodes):
    """Check if subset is a proper subset of {1, ..., num_remote_nodes}."""
    if not isinstance(subset, set):
        raise TypeError("subset must be a set.")
    if not subset.issubset(set(range(1, num_remote_nodes + 1))):
        raise ValueError("subset must be a subset of {1, 2, ..., num_remote_nodes}.")