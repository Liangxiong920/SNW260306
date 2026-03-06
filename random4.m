% =========================================================================
% NUMERICAL INVESTIGATION OF RANDOM MIXED STATES AND WITNESS BOUNDS
%
% This script numerically investigates different mixed states by generating 
% 50 independent samples using MATLAB. Each random density matrix is 
% constructed as a finite convex combination of independently sampled 
% Haar-distributed pure states. Specifically, for each sample, 10000 independent 
% Haar-random pure states |\psi_j> are generated via QR decomposition of 
% complex Ginibre matrices. Random positive weights are assigned and the 
% resulting operator is normalized to unit trace, ensuring a valid density 
% matrix. The optimal witness coefficient \lambda_{4} and upper bounds 
% of the mixed state are derived and analyzed.
% =========================================================================

% Set random seed for reproducibility
rng(1);

% Define system dimensions
m = 3;
n = 3;

% Define the number of density matrices to generate
num_density_matrices = 50;

% Initialize arrays to store function values for each density matrix
lambda_4_values = zeros(num_density_matrices, 1);
theta_4_values = zeros(num_density_matrices, 1);
%zeta_4_values = zeros(num_density_matrices, 1);
eta_4_values = zeros(num_density_matrices, 1);
P_4_values = zeros(num_density_matrices, 1);

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
    s5=sorted_singular_values(5); s6=sorted_singular_values(6); s7=sorted_singular_values(7); s8=sorted_singular_values(8);
    s9=sorted_singular_values(9); %s10=sorted_singular_values(10); s11=sorted_singular_values(11); s12=sorted_singular_values(12);
    %s13=sorted_singular_values(13); s14=sorted_singular_values(14); s15=sorted_singular_values(15); s16=sorted_singular_values(16);
%     s17=sorted_singular_values(17); s18=sorted_singular_values(18); s19=sorted_singular_values(19); s20=sorted_singular_values(20);
%     s21=sorted_singular_values(21); s22=sorted_singular_values(22); s23=sorted_singular_values(23); s24=sorted_singular_values(24);
%     s25=sorted_singular_values(25);

    % Matrix M3
    M3= [s1 s2 s4; s3 s6 s7; s5 s8 s9];

    % Calculate eigenvalues
    eigenvalues = eig(M3);

    % Find the maximum eigenvalue
    max_eigenvalue = max(eigenvalues);
    lambda_4 = max_eigenvalue;
    theta_4 = s1+s2+s4-s9*(1-(2*(s5+s8)/(s1+s2+s4-2*s9+sqrt((s1+s2+s4)^2-4*s9*((s1+s2+s4)-(s5+s8+s9))))));
    %zeta_4 = s1+s2+s4-s9*(1-sqrt((s5+s8)/(s1+s2+s4-s9)));
    eta_4 = s1+s2+s4-s9*(1-sqrt((s5+s8+s9)/(s1+s2+s4)));
    P_4 = s1+s2+s4;

    % Store function values
    lambda_4_values(k) = lambda_4;
    theta_4_values(k) = theta_4;
    %zeta_4_values(k) = zeta_4;
    eta_4_values(k) = eta_4;
    P_4_values(k) = P_4;
end


figure;
plot(1:num_density_matrices, lambda_4_values, 'o-', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '\lambda_4', 'Color', [0 0.4470 0.7410]);
hold on;
plot(1:num_density_matrices, theta_4_values, 'x--', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '\theta_4', 'Color', [0.8500 0.3250 0.0980]);
%plot(1:num_density_matrices, zeta_4_values, 's-.', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '\zeta_4', 'Color', [0.9290 0.6940 0.1250]);
plot(1:num_density_matrices, eta_4_values, 'd:', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '\eta_4', 'Color', [0.4940 0.1840 0.5560]);
plot(1:num_density_matrices, P_4_values, '*-', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'P_4', 'Color', [0.4660 0.6740 0.1880]);
xlabel('Density Matrix Index');
ylabel('Coefficient and Bounds Values');
title('The optimal witnesses coefficient \lambda_{4} and upper bounds');
legend show;
grid on;

