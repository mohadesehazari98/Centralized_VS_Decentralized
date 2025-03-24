function [Rate_out] = Rate_2D(q_BSM, q_Fuse, N, delta_t, L_0_in, m, k_max, t)
% Ref: Mohadeseh Azari / Department of Informatics and Networked Systems / 
% School of Computing and Information / University of Pittsburgh / 
% Pittsburgh,PA / moa125@pitt.edu
% This function computes the average rate for distributing N-qubit GHZ entanglement 
% through an m-level 2D repeater approach, utilizing an either centralized or 
% decentralized switch for the generation of the parent entanglement.

% Parameters:
% q_BSM: The probability of successful Bell state measurement for each qubit.
% N: Number of qubits in the entangled state.
% delta_t: The time step of the protocol (in seconds).
% L_0_in: The final distance between neighboring nodes to which the 
% entangled state is teleported (in kilometers).
% m: The number of 2D repeater generation (children) in the network.

q_link_values = link_gen_prob(t, L_0_in, m, N);

% Initialize an array to store the expected time for each link (ETmax)
ETmax = zeros(size(q_link_values));

for idx = 1:length(q_link_values)
    
    % Get the current link value (q_link) for this iteration
    q_link = q_link_values(idx);
    
    % Initialize a variable to accumulate the sum for the expected time (ETmax)
    sum_ETmax = 0;
    
    for k = 1:k_max

        % The probability that the parent hasn't finished by time k*delta_t
        % is 1 - F_Tmax_k.
        term = 1 - F_T_max(k, m, N, q_BSM, q_Fuse, q_link, t);
        
        % Accumulate the term to the total sum of expected times
        sum_ETmax = sum_ETmax + term;
    end
    
    % Store the result for the expected time
    ETmax(idx) = delta_t * sum_ETmax;
end

% Calculate the rate of entanglement distribution based on the expected
% time to have all the parents ready and having all final BSM succeed
% simultaneously 
Rate_out = q_BSM^(N*(N-1)/2) ./ ETmax;

end



