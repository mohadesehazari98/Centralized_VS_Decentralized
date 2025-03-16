function [Rate_out] = Rate_Cent(q_BSM, N, delta_t, L_0_in)
% This function calculates the average rate of distributing N-qubit GHZ
% entangled state using a centralized switch.

% Parameters:
% q_BSM: The Probability of successful Bell state measurements.
% N: The number of qubits involved in the distributed entangled state.
% delta_t: The time step for each calculation (usually in seconds).
% L_0_in: The final distance between neighboring nodes to which the 
% entangled state is teleported (in kilometers).

% Step 1: Adjust the distance to the central switch based on geometric
% scaling factors.
L_0 = L_0_in ./ (2 * sin(pi / N));

% Step 2: Define constants for link efficiency and attenuation.
etha_c = 0.95;  % coupling efficiency (emission of the photon from the memory qubit)
L_att = 20;  % Link attenuation (distance over which signal weakens by 1/e, in km)

% Step 3: Calculate the link probability (q_link) based on the initial link length and attenuation
q_link = 0.5 .* etha_c^2 .* exp(- L_0 ./ L_att); 
% This formula calculates the probability of successful transmission over a distance.

ET = zeros(size(q_link)); 

% Step 4: Compute the expected time (E[T]) for each link value (q_link)
for idx = 1:length(q_link)
    
    % Get the current q_link value for this iteration
    q_link_values = q_link(idx);
    
    sum_term = 0;
    
    % Step 5: Compute the summation term over the number of qubits (N)
    for j = 1:N
        % Compute the binomial coefficient (N choose j) - represents the number of ways
        % to choose j successes out of N possible outcomes.
        binomial_coeff = nchoosek(N, j); 
        
        % Compute the fraction term based on the link probability (q_link)
        fraction_term = 1 / (1 - (1 - q_link_values)^j); 
        
        % Add or subtract this term to the running total (sum_term)
        sum_term = sum_term + (-1)^(j+1) * binomial_coeff * fraction_term;
    end
    
    % Step 6: Compute the expected time (E[T]) using the formula
    ET(idx) = (delta_t / q_BSM.^N) * sum_term; 
end

% Step 7: Calculate the overall rate by taking the inverse of the expected time (E[T])
Rate_out = 1 ./ ET; 

end
