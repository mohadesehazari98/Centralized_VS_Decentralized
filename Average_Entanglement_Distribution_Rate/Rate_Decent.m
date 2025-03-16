function Rate_out = Rate_Decent(q_BSM, q_Fuse, N, delta_t, L_0_in, k_max)
% Rate_Decent calculates the average GHZ entanglement distribution rate
% for an N-qubit system using a decentralized switching approach.
%
% INPUTS:
%   q_BSM    - Success probability of a Bell-state measurement (BSM)
%   q_Fuse   - Success probability of a fusion operation
%   N        - Number of qubits/nodes involved in GHZ entanglement
%   delta_t  - Time duration of a single attempt (time step in s)
%   L_0_in   - The final distance between neighboring nodes to which the 
%              entangled state is teleported (in km)
%
% OUTPUT:
%   Rate     - Average entanglement distribution rate (inverse of E[T_max])

%----------------------------
% Step 1: Adjust Link Length
%----------------------------
% Adjust L_0 based on geometric scaling.
L_0 = L_0_in ./ 2;

%-------------------------------------
% Step 2: Define Link Efficiency Model
%-------------------------------------
etha_c = 0.98;      % Coupling efficiency (e.g., probability of emitting photon from memory qubit)
L_att  = 20;        % Attenuation length (distance over which signal weakens by 1/e, in km)

% Compute overall link success probability (includes coupling loss and fiber loss)
% q_link = 0.5 * eta_c^2 * exp(-L_0 / L_att);
q_link_values = 0.5 .* etha_c^2 .* exp(-L_0 ./ L_att);

% Initialize an array to store the expected time (E[T]) for each link probability (q_link)
E_T = zeros(size(q_link_values)); 

for idx = 1:length(q_link_values)
    
    % Get the current q_link value for this iteration
    q_link = q_link_values(idx);

    % Initialize arrays to store CDF values
    F_n_i = zeros(k_max,1);  % CDF of a single n_i (individual system completion time)
    F_n   = zeros(k_max,1);  % CDF of the maximum completion time: n = max{n_1, ..., n_N}
    
    %---------------------------------------------------------
    % Step 3: Compute CDF of individual n_i using convolution
    %---------------------------------------------------------
    % Each n_i = n_BSM * n_dist, where:
    % - n_BSM is geometric(q_BSM)
    % - n_dist has CDF: [1 - (1 - q_link)^v]^2
    
    for k = 1:k_max
        sum_u = 0; 
        for u = 1:k
            % Probability that n_BSM = u (geometric distribution)
            p_nBSM = (1 - q_BSM)^(u - 1) * q_BSM;
    
            % For fixed n_BSM = u, we want P(n_dist <= floor(k/u))
            L = floor(k / u);  % Maximum possible n_dist
            p_nDist_CDF = (1 - (1 - q_link)^L)^2;  % CDF of n_dist at L
    
            % Contribution to overall CDF of n_i
            sum_u = sum_u + p_nBSM * p_nDist_CDF;
        end
        F_n_i(k) = sum_u;              % CDF of n_i up to time slot k
        F_n(k) = F_n_i(k)^N;           % CDF of max{n_i} using independence: F_n_i(k)^N
    end
    
    %-------------------------------------------------
    % Step 4: Compute E[n] (Expected max completion time)
    %-------------------------------------------------
    E_n = 0;
    for k = 1:k_max
        E_n = E_n + 1 - F_n(k);
    end
    
    %----------------------------------------
    % Step 5: Compute Expected Total Time E[T]
    %----------------------------------------
    % Multiply by time step delta_t and account for fusion operations
    % Each fusion step must succeed independently for all N fusions
    E_T(idx) = (E_n .* delta_t) ./ (q_Fuse)^N;
end
%-------------------------------------
% Step 6: Compute Final Entanglement Rate
%-------------------------------------
Rate_out = 1 ./ E_T;  % Inverse of average time gives the entanglement rate

end