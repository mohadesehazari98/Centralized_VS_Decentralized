function result = F_T_max(k, m, N, q_BSM, q_link)
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
