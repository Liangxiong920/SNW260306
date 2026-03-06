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
% matrix. The optimal witness coefficient \lambda_{5} and upper bounds 
% of the mixed state are derived and analyzed.
% =========================================================================

% Set random seed for reproducibility
rng(1);

% Define system dimensions
m = 4;
n = 4;

% Define the number of density matrices to generate
num_density_matrices = 50;

% Initialize arrays to store function values for each density matrix
lambda_5_values = zeros(num_density_matrices, 1);
theta_5_values = zeros(num_density_matrices, 1);
%zeta_5_values = zeros(num_density_matrices, 1);
eta_5_values = zeros(num_density_matrices, 1);
P_5_values = zeros(num_density_matrices, 1);

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
    s9=sorted_singular_values(9); s10=sorted_singular_values(10); s11=sorted_singular_values(11); s12=sorted_singular_values(12);
    s13=sorted_singular_values(13); s14=sorted_singular_values(14); s15=sorted_singular_values(15); s16=sorted_singular_values(16);
%     s17=sorted_singular_values(17); s18=sorted_singular_values(18); s19=sorted_singular_values(19); s20=sorted_singular_values(20);
%     s21=sorted_singular_values(21); s22=sorted_singular_values(22); s23=sorted_singular_values(23); s24=sorted_singular_values(24);
%     s25=sorted_singular_values(25);

    % Matrix M4
    M4= [s1 s2 s4 s6; s3 s8 s9 s11; s5 s10 s13 s14; s7 s12 s15 s16];

    % Calculate eigenvalues
    eigenvalues = eig(M4);

    % Find the maximum eigenvalue
    max_eigenvalue = max(eigenvalues);
    lambda_5 = max_eigenvalue;
    %iota_5=(s1^2+s2^2+s3^2+s4^2+s5^2+s6^2+s7^2+s8^2+s9^2+(s10)^2+(s11)^2+(s12)^2+(s13)^2+(s14)^2+(s15)^2+(s16)^2)^(1/2);
    theta_5 = s1+s2+s4+s6-s16*(1-(2*(s7+s12+s15)/(s1+s2+s4+s6-2*s16+sqrt((s1+s2+s4+s6)^2-4*s16*((s1+s2+s4+s6)-(s7+s12+s15+s16))))));
    %zeta_5 = s1+s2+s4+s6-s16*(1-sqrt((s7+s12+s15)/(s1+s2+s4+s6-s16)));
    eta_5 = s1+s2+s4+s6-s16*(1-sqrt((s7+s12+s15+s16)/(s1+s2+s4+s6)));
    P_5 = s1+s2+s4+s6;

    % Store function values
    lambda_5_values(k) = lambda_5;
    %iota_5_values(k) = iota_5;
    theta_5_values(k) = theta_5;
    %zeta_5_values(k) = zeta_5;
    eta_5_values(k) = eta_5;
    P_5_values(k) = P_5;
end

% Plot function values
% figure;
% plot(1:num_density_matrices, lambda_6_values, 'o-', 'DisplayName', '\lambda_6');
% hold on;
% plot(1:num_density_matrices, theta_6_values, 'x-', 'DisplayName', '\theta_6');
% plot(1:num_density_matrices, zeta_6_values, 's-', 'DisplayName', '\zeta_6');
% plot(1:num_density_matrices, eta_6_values, 'd-', 'DisplayName', '\eta_6');
% plot(1:num_density_matrices, P_6_values, '*-', 'DisplayName', 'P_6');
% xlabel('Density Matrix Index');
% ylabel('Function Value');
% title('Function Values for Random Density Matrices');
% legend show;
% grid on;

figure;
plot(1:num_density_matrices, lambda_5_values, 'o-', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', '\lambda_5', 'Color', [0 0.4470 0.7410]);
hold on;
plot(1:num_density_matrices, theta_5_values, 'x--', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', '\theta_5', 'Color', [0.8500 0.3250 0.0980]);
%plot(1:num_density_matrices, zeta_5_values, 's-.', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', '\zeta_5', 'Color', [0.9290 0.6940 0.1250]);
plot(1:num_density_matrices, eta_5_values, 'd:', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', '\eta_5', 'Color', [0.4940 0.1840 0.5560]);
plot(1:num_density_matrices, P_5_values, '*-', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'P_5', 'Color', [0.4660 0.6740 0.1880]);
% Plot iota_5
plot(1:num_density_matrices, iota_5_values, 'v-.', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', '\iota_5', 'Color', [0.3010 0.7450 0.9330]);
xlabel('Density Matrix Index');
ylabel('Coefficient and Bounds Values');
title('The optimal witnesses coefficient \lambda_{5} and upper bounds');
legend show;
grid on;

