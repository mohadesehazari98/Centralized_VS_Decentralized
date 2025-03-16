function [Rate_out] = Rate_2D_Cent(q_BSM, N, delta_t, L_0_in, m, k_max)
% This function computes the average rate for distributing N-qubit GHZ entanglement 
% through an m-level 2D repeater approach, utilizing a centralized switch for 
% the generation of the parent entanglement.

% Parameters:
% q_BSM: The probability of successful Bell state measurement for each qubit.
% N: Number of qubits in the entangled state.
% delta_t: The time step of the protocol (in seconds).
% L_0_in: The final distance between neighboring nodes to which the 
% entangled state is teleported (in kilometers).
% m: The number of 2D repeater generation (children) in the network.

% Adjust the initial link length L_0 based on the distance scaling factor
L_0 = L_0_in ./ (2^m * 2 * sin(pi / N));

% Define constants for attenuation and link efficiency
etha_c = 0.95;  % coupling efficiency (emission of the photon from the memory qubit)
L_att = 20;  % Link attenuation in kilometers (distance over which signal weakens by 1/e)

% Calculate the link values based on the initial link length and attenuation
q_link_values = 0.5 .* etha_c^2 .* exp(- L_0 ./ L_att);

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
        term = 1 - F_T_max(k,m,N,q_BSM,q_link);
        
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



