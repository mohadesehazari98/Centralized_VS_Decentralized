function Rate = Rate_2D_Decent(q_BSM, q_Fuse, N, delta_t, L_0_in, m, k_max)
% Rate_2D_Decent computes the average entanglement distribution rate for 
% an N-qubit GHZ state in a 2D repeater architecture with decentralized fusion.
%
% INPUTS:
%   q_BSM    - Success probability of Bell-state measurement
%   q_Fuse   - Success probability of fusion operation
%   N        - Number of qubits involved in GHZ entanglement
%   delta_t  - Time step duration for one attempt
%   L_0_in   - Initial physical link length between neighboring nodes (in km)
%   m        - Hierarchical level of fusion (depth of 2D construction)
%
% OUTPUT:
%   Rate     - Entanglement distribution rate (1 / E[T_max])

%-----------------------------------------
% Step 1: Setup summation limits and tolerances
%-----------------------------------------

%------------------------------------------------------
% Step 2: Adjust link length based on hierarchical level
%------------------------------------------------------
% Each hierarchical level reduces link length by a factor of 2^m
% sqrt(3)/3 is a geometric scaling factor for 2D triangular grids
L_0 = (1 / (2^m)) * L_0_in;

%------------------------------------------------
% Step 3: Link attenuation and link success model
%------------------------------------------------
etha_c = 0.95;        % Coupling efficiency
L_att  = 20;          % Attenuation length of fiber (in km)

% Calculate link success probability q_link
% Includes coupling loss and fiber attenuation
q_link = 0.5 * etha_c^2 .* exp(-L_0 / L_att);

%---------------------------------------
% Step 4: Initialize expected time array
%---------------------------------------
% E_Tmax stores the expected max completion time for each link setting
E_Tmax = zeros(size(q_link)); 

%--------------------------------------------------------
% Step 5: Main loop to compute E[T_max] using tail-sum
%--------------------------------------------------------
for idx = 1:length(q_link)
    q_link_values = q_link(idx);   % Pick current link value

    % Initialize E[T_max] for current link
    for k = 1:k_max
        FT = 0;  % CDF F_{T_max}(k)

        % Outer sum over u: fusion success geometric distribution
        for u = 1:k
            % Probability that fusion success takes u rounds
            P_mF = (1 - q_Fuse^N)^(u - 1) * (q_Fuse^N);

            % Maximum index for BSM and distribution attempt
            k_inner = floor(k / u);
            F_n_i = 0;  % CDF of individual system n_i

            % Inner sum over v: BSM success geometric distribution
            for v = 1:k_inner
                % Probability that BSM takes v rounds
                P_BSM = (1 - q_BSM)^(v - 1) * q_BSM;

                % Probability that distribution succeeds within ⌊k_inner / v⌋ rounds
                dist_term = 1 - (1 - q_link_values)^floor(k_inner / v);
                P_dist = dist_term^2;

                % CDF of individual n_i
                F_n_i = F_n_i + P_BSM * P_dist;
            end

            % Contribution to F_T for current u
            FT = FT + P_mF * F_n_i^N;
        end

        % Total CDF of T_max at time k
        F_Tmax = FT^(N^m);

        % Tail-sum: Add P(T_max ≥ k) = 1 - F_Tmax
        E_Tmax(idx) = E_Tmax(idx) + (1 - F_Tmax);
    end

    % Multiply by time step to get E[T_max] in time units
    E_Tmax(idx) = delta_t * E_Tmax(idx);
end

%---------------------------------------------------
% Step 6: Compute final entanglement distribution rate
%---------------------------------------------------
% Rate is normalized by BSM success across all links: q_BSM^(N(N-1)/2)
Rate = q_BSM^(N * (N - 1) / 2) ./ E_Tmax;

end
