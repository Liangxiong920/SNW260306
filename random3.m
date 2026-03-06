
% =========================================================================
% NUMERICAL INVESTIGATION OF RANDOM MIXED STATES AND WITNESS BOUNDS
%
% This script numerically investigates mixed states by generating 50 
% independent random density matrices. 
%
% Key details of the simulation:
% - State Generation: Each density matrix is constructed as a finite convex 
%   combination of 10,000 independent Haar-distributed pure states.
% - Methodology: The pure states are generated via the QR decomposition of 
%   complex Ginibre matrices. They are assigned random positive weights and 
%   normalized to unit trace to ensure valid density matrices.
% - Scope: The analysis focuses on bipartite systems with Schmidt numbers 
%   not exceeding k=2.
% - Objective: For each mixed state, the script calculates the optimal 
%   witness coefficient (\lambda_3) and its corresponding upper bounds, 
%   plotting the results for comparison.
% =========================================================================

% Set random seed for reproducibility
rng(1);

% Define system dimensions
m = 2;
n = 2;

% Define the number of density matrices to generate
num_density_matrices = 50;

% Initialize arrays to store function values for each density matrix
lambda_3_values = zeros(num_density_matrices, 1);
theta_3_values = zeros(num_density_matrices, 1);
%zeta_3_values = zeros(num_density_matrices, 1);
eta_3_values = zeros(num_density_matrices, 1);
P_3_values = zeros(num_density_matrices, 1);

for k = 1:num_density_matrices
    % Initialize the density matrix
    rho = zeros(m*n, m*n);

    % Define the number of pure states to generate
    num_pure_states = 10000;

    for j = 1:num_pure_states
        % Generate a random complex matrix
        A = randn(m*n) + 1i*randn(m*n);

        % Perform QR decomposition on the matrix to obtain a random unitary matrix
        [Q, ~] = qr(A);

        % Select the first column vector of the unitary matrix as the pure state
        psi = Q(:, 1);

        % Generate random probability weight
        p = rand();

        % Accumulate the weighted pure state density matrix
        rho = rho + p * (psi * psi');
    end

    % Normalize the density matrix
    rho = rho / trace(rho);

    % Expand the density matrix rho row by row into a row vector
    rho_vector = rho(:).';

    % Initialize the correlation matrix
    correlation_matrix = zeros(size(rho));

    % Calculate each element of the correlation matrix
    for i = 1:size(rho,1)
        for j = 1:size(rho,2)
            unit_vector_i = zeros(1, size(rho,1));
            unit_vector_i(i) = 1; % First subsystem unit orthogonal basis column vector
            unit_vector_j = zeros(1, size(rho,2));
            unit_vector_j(j) = 1; % Second subsystem unit orthogonal basis column vector

            % Calculate the tensor product result
            tensor_product = kron(unit_vector_i, unit_vector_j);

            % Calculate the elements of the correlation matrix
            correlation_matrix(i,j) = rho_vector * tensor_product(:);
        end
    end

    % Calculate singular values
    singular_values = svd(correlation_matrix);

    % Sort singular values in descending order
    sorted_singular_values = sort(singular_values, 'descend');
    s1=sorted_singular_values(1); s2=sorted_singular_values(2); s3=sorted_singular_values(3); s4=sorted_singular_values(4);
    


    % Matrix M
    M2= [s1 s2; s3 s4];

    % Calculate eigenvalues
    eigenvalues = eig(M2);

    % Find the maximum eigenvalue
    max_eigenvalue = max(eigenvalues);
    lambda_3 = max_eigenvalue;
    theta_3 = s1+s2-s4*(1-(2*(s3)/(s1+s2-2*s4+sqrt((s1+s2)^2-4*s4*((s1+s2)-(s3+s4))))));
    %zeta_3 = s1+s2-s4*(1-sqrt((s3)/(s1+s2-s4)));
    eta_3 = s1+s2-s4*(1-sqrt((s3+s4)/(s1+s2)));
    P_3 = s1+s2;

    % Store function values
    lambda_3_values(k) = lambda_3;
    theta_3_values(k) = theta_3;
    %zeta_3_values(k) = zeta_3;
    eta_3_values(k) = eta_3;
    P_3_values(k) = P_3;
end


figure;
plot(1:num_density_matrices, lambda_3_values, 'o-', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', '\lambda_3', 'Color', [0 0.4470 0.7410]);
hold on;
plot(1:num_density_matrices, theta_3_values, 'x--', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', '\theta_3', 'Color', [0.8500 0.3250 0.0980]);
%plot(1:num_density_matrices, zeta_3_values, 's-.', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', '\zeta_3', 'Color', [0.9290 0.6940 0.1250]);
plot(1:num_density_matrices, eta_3_values, 'd:', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', '\eta_3', 'Color', [0.4940 0.1840 0.5560]);
plot(1:num_density_matrices, P_3_values, '*-', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', 'P_3', 'Color', [0.4660 0.6740 0.1880]);
xlabel('Density Matrix Index');
ylabel('Coefficient and Bounds Values');
title('The optimal witnesses coefficient \lambda_{3} and upper bounds');
legend show;
grid on;

