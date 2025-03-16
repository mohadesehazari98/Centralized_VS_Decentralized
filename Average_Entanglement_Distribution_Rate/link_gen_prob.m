function q_link_values = link_gen_prob(t, L_0_in, m, N)
% link_gen_prob computes the probability of successful link generation q_link
% based on repeater type (centralized/decentralized), system parameters, and geometry.
%
% INPUTS:
%   t        - 'Centralized' or 'Decentralized' (string or char)
%   L_0_in   - Base distance (before scaling), in km
%   m        - Hierarchical level in the repeater structure
%   N        - Number of nodes/qubits
%
% OUTPUT:
%   q_link_values - Probability of successful link generation (array)

% Use strcmpi for case-insensitive string comparison
if strcmpi(t, 'Centralized')
    % For Centralized GHZ distribution, effective link length:
    % L = L_0_in / (2^m * 2 * sin(pi / N))
    L_0 = L_0_in ./ (2^m * 2 * sin(pi / N));
    
elseif strcmpi(t, 'Decentralized')
    % For Decentralized GHZ distribution, effective link length:
    % L = L_0_in / (2^m * 2)
    L_0 = L_0_in ./ (2^m * 2);
    
else
    error('Invalid type t. Use either ''Centralized'' or ''Decentralized''.');
end

% Define constants for attenuation and coupling efficiency
etha_c = 0.95;     % Coupling efficiency (photon emission from memory qubit)
L_att = 20;        % Attenuation length (km) — signal weakens by 1/e over this length

% Compute q_link based on attenuation and distance
q_link_values = 0.5 .* etha_c^2 .* exp(- L_0 ./ L_att);

end
