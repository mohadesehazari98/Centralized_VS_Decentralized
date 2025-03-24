clear; clc;
% This script plots the average rate of multipartite entanglement distribution
% as a function of end-node distance in both centralized and decentralized schemes.
% The GHZ state involves N qubits/nodes.
% Parameters
q_BSM_values = 0.98;                 % Bell-state measurement success probability
q_Fuse_values = q_BSM_values;        % Fusion success probability (here set equal to q_BSM)
N = [3, 4, 5];                       % Number of qubits/nodes in GHZ entanglement
delta_t = 1;                         % Time step duration
L_lin = linspace(0.001, 1, 100).^2;  % Square a linearly spaced vector to push points toward zero
L_0_in = L_lin * 150;                % Scale to desired range [0, 200]
k_max = 2000;                        % Truncation parameter for infinite summation
% Plotting
figure(1);
hold on;
color_order = get(gca, 'ColorOrder');
num_colors = size(color_order, 1);
for i = 1:length(N)
    % Compute rates
    rate_Cent = Rate_Cent(q_BSM_values, N(i), delta_t, L_0_in);
    rate_Decent = Rate_Decent(q_BSM_values, q_Fuse_values, N(i), delta_t, L_0_in, k_max);
    % Select color for both plots
    color_idx = mod(i-1, num_colors) + 1;
    current_color = color_order(color_idx, :);
    % Plot centralized scheme
    plot(L_0_in, rate_Cent, '-', 'DisplayName', ['Centralized, N = ' num2str(N(i))], ...
        'Color', current_color, 'LineWidth', 2);
    % Plot decentralized scheme with same color, dashed line
    plot(L_0_in, rate_Decent, '--', 'DisplayName', ['Decentralized, N = ' num2str(N(i))], ...
        'Color', current_color, 'LineWidth', 2);
end
% Final touches
legend show;
xlabel('End-node Distance (km)', 'FontSize', 12);
ylabel('Entanglement Distribution Rate', 'FontSize', 12);
title(['Entanglement Distribution: Centralized vs. Decentralized (q_{BSM} = ' num2str(q_BSM_values) ')'], 'FontSize', 14);
grid on;

%% 
% This script plots the average rate of multipartite entanglement distribution
% as a function of end-node distance in both centralized and decentralized schemes.
% The N is constant at 3 and q_BSM is changing.
% Parameters
q_BSM_values = [0.5, 0.7, 0.9, 1];   % Bell-state measurement success probability
q_Fuse_values = q_BSM_values;        % Fusion success probability (here set equal to q_BSM)
N = 4;                               % Number of qubits/nodes in GHZ entanglement
delta_t = 1;                         % Time step duration
L_lin = linspace(0.001, 1, 100).^2;  % Square a linearly spaced vector to push points toward zero
L_0_in = L_lin * 150;                % Scale to desired range [0, 200]
k_max = 2000;                        % Truncation parameter for infinite summation
% Plotting
figure(2);
hold on;
color_order = get(gca, 'ColorOrder');
num_colors = size(color_order, 1);
for i = 1:length(q_BSM_values)
    % Compute rates
    rate_Cent = Rate_Cent(q_BSM_values(i), N, delta_t, L_0_in);
    rate_Decent = Rate_Decent(q_BSM_values(i), q_Fuse_values(i), N, delta_t, L_0_in, k_max);
    % Select color for both plots
    color_idx = mod(i-1, num_colors) + 1;
    current_color = color_order(color_idx, :);
    % Plot centralized scheme
    plot(L_0_in, rate_Cent, '-', 'DisplayName', ['Centralized, q_{BSM} = ' num2str(q_BSM_values(i))], ...
        'Color', current_color, 'LineWidth', 2);
    % Plot decentralized scheme with same color, dashed line
    plot(L_0_in, rate_Decent, '--', 'DisplayName', ['Decentralized, q_{BSM} = ' num2str(q_BSM_values(i))], ...
        'Color', current_color, 'LineWidth', 2);
end
% Final touches
legend show;
xlabel('End-node Distance (km)', 'FontSize', 12);
ylabel('Entanglement Distribution Rate', 'FontSize', 12);
title(['Entanglement Distribution: Centralized vs. Decentralized (N = ' num2str(N) ')'], 'FontSize', 14);
grid on;

%%
% This script plots the average rate of multipartite entanglement distribution
% in a 2D scheme as a function of end-node distance in both centralized and 
% decentralized schemes.
% Parameters
q_BSM_values = 0.98;                 % Bell-state measurement success probability
q_Fuse_values = q_BSM_values;        % Fusion success probability (here set equal to q_BSM)
N = [3, 4, 5];                       % Number of qubits/nodes in GHZ entanglement
delta_t = 1;                         % Time step duration
L_lin = linspace(0.001, 1, 100).^2;  % Push points toward zero
L_0_in = L_lin * 150;                % Scale to desired range [0, 200]
k_max = 2000;                        % Truncation parameter
m = [1, 2, 3];                       % Nesting levels
% Plotting
figure(3);
color_order = get(gca, 'ColorOrder');
num_colors = size(color_order, 1);
for i = 1:length(m)
    subplot(1, length(m), i);
    cla; hold on;
    for j = 1:length(N)
        % Compute rates
        rate_Cent = Rate_2D(q_BSM_values, q_Fuse_values, N(j), delta_t, L_0_in, m(i), k_max, 'Centralized');
        rate_Decent = Rate_2D(q_BSM_values, q_Fuse_values, N(j), delta_t, L_0_in, m(i), k_max, 'Decentralized');
        % Select color
        color_idx = mod(j-1, num_colors) + 1;
        current_color = color_order(color_idx, :);
        % Plot centralized
        plot(L_0_in, rate_Cent, '-', 'DisplayName', ['Centralized, N = ' num2str(N(j))], ...
            'Color', current_color, 'LineWidth', 2);
        % Plot decentralized
        plot(L_0_in, rate_Decent, '--', 'DisplayName', ['Decentralized, N = ' num2str(N(j))], ...
            'Color', current_color, 'LineWidth', 2);
    end
    legend show;
    xlabel('End-node Distance (km)', 'FontSize', 12);
    ylabel('Entanglement Distribution Rate', 'FontSize', 12);
    title(['Nesting Level m = ' num2str(m(i))], 'FontSize', 14);
    grid on;
end
