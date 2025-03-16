function result = F_T_max(k, m, N, q_BSM, q_Fuse, q_link, t)
% F_T_max recursively computes the CDF of the maximum completion time T_max
% for a hierarchical (multi-level) quantum repeater system.
%
% INPUTS:
%   k      - integer value at which F_T_max(k) is evaluated
%   m      - recursion level (e.g., m=1 means one level above physical links)
%   N      - number of qubits (nodes) at each level
%   q_BSM  - success probability of Bell State Measurement (BSM)
%   q_link - success probability of stablishing elementary link
%
% OUTPUT:
%   result - the computed value of F_T_max(k)

% Base case: when m = 1, we compute the first-level CDF using the link-level statistics
if m == 1
    if strcmpi(t, 'Centralized')

        F_T_k = 0;  % Initialize the sum for CDF of T_max at level m=1
        for u = 1:k
            % Probability that state teleportation succeeds at the u-th attempt:
            % P(n_teleport = u) = (1 - q_BSM^N)^(u-1) * q_BSM^N
            p_n1 = (1 - q_BSM^N)^(u-1) * q_BSM^N;
            
            % For each u, determine how many link-level attempts are required
            exponent = floor(k/u);
            
            % CDF of max of N link-level attempts: (F_n_i)^N
            % where F_n_i = P(link success ≤ exponent) = [1 - (1 - q_link)^exponent]
            F_n2_max = (1 - (1 - q_link)^exponent)^N;
            
            % Add weighted contribution to total CDF
            F_T_k = F_T_k + p_n1 * F_n2_max;
        end
        % Final CDF is raised to power N due to max over N identical systems
        result = F_T_k^N;

    elseif strcmpi(t, 'Decentralized')

        F_T_k = 0;  % CDF F_{T_max}(k)
        for u = 1:k
            % Probability that fusion success takes u rounds
            P_mF = (1 - q_Fuse^N)^(u - 1) * (q_Fuse^N);

            % Maximum index for BSM and distribution attempt
            exponent = floor(k / u);
            F_n_i = 0;  % CDF of individual system n_i

            % Inner sum over v
            for v = 1:exponent
                % Probability that BSM takes v rounds
                P_BSM = (1 - q_BSM)^(v - 1) * q_BSM;

                % Probability that distribution succeeds within ⌊exponent / v⌋ rounds
                P_dist = (1 - (1 - q_link)^floor(exponent / v))^2;

                % CDF of individual n_i
                F_n_i = F_n_i + P_BSM * P_dist;
            end

            % Contribution to F_T_k for current u
            F_T_k = F_T_k + P_mF * F_n_i^N;
        end

        % Total CDF of T_max at time k
        result = F_T_k^N;
        
    else
        error('Invalid type t. Use either ''Centralized'' or ''Decentralized''.');
    end
else
    % Recursive case: m > 1, compute based on lower level F_T_max
    F_T_k = 0;  % Initialize sum for current level
    
    % Loop over number of fusion attempts u = 1 to k
    for u = 1:k
        % Probability that fusion of parent GHZ succeeds at u-th attempt
        % In general, q_fuse = q_BSM^(N choose 2) = q_BSM^(N*(N-1)/2)
        p_n1 = (1 - q_BSM^(N*(N-1)/2))^(u-1) * q_BSM^(N*(N-1)/2);
        
        % Compute effective trials needed at lower level
        exponent = floor(k/u);
        
        % Recursive call to compute lower-level F_T_max
        F_T_k = F_T_k + p_n1 * F_T_max(exponent, m-1, N, q_BSM, q_link);
    end
    
    % Raise the sum to the N-th power due to maximum over N systems
    result = F_T_k^N;
end
end
