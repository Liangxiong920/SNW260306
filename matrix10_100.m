% =========================================================================
% NUMERICAL INVESTIGATION OF RANDOM MIXED STATES AND WITNESS BOUNDS
%
% This script numerically investigates different mixed states by generating 
% 50 independent samples using MATLAB. Each random density matrix is 
% constructed as a finite convex combination of independently sampled 
% Haar-distributed pure states. Specifically, for each sample, 100 independent 
% Haar-random pure states |\psi_j> are generated via QR decomposition of 
% complex Ginibre matrices. Random positive weights are assigned and the 
% resulting operator is normalized to unit trace, ensuring a valid density 
% matrix. The optimal witness coefficient \lambda_{11} and upper bounds 
% of the mixed state are derived and analyzed.
% =========================================================================

% Set random seed for reproducibility
rng(1);

% Define system dimensions
m = 10;
n = 10;

% Define the number of density matrices to generate
num_density_matrices = 50;

% Initialize arrays to store function values for each density matrix
lambda_11_values = zeros(num_density_matrices, 1);
theta_11_values = zeros(num_density_matrices, 1);
eta_11_values = zeros(num_density_matrices, 1);
P_11_values = zeros(num_density_matrices, 1);

for k = 1:num_density_matrices
    % Initialize the density matrix
    rho = zeros(m*n, m*n);

    % Define the number of pure states to generate
    num_pure_states = 100;

     for j = 1:num_pure_states
        % Generate a random complex matrix
        A = randn(m*n) + 1j*randn(m*n);

        % Perform QR decomposition on the matrix to obtain a random unitary matrix
        [Q, ~] = qr(A);

        % Select the first column vector of the unitary matrix as the pure state
        psi = Q(:, 1);

        % Generate random probability weight
        p = rand();

        % Accumulate the weighted pure state density matrix
        rho = rho + p * (psi * ctranspose(psi));
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

    % Calculate singular value decomposition (SVD)
    [U, S, V] = svd(correlation_matrix);

    % Extract singular values
    singular_values = diag(S);

    % Sort singular values in descending order
    [singular_values_sorted, idx] = sort(singular_values, 'descend');

    % Assign the first 100 singular values to s1, s2, ..., s100
    for i = 1:100
    eval(['s' num2str(i) ' = singular_values_sorted(i);']);
    end


%iato_values = zeros(num_density_matrices, 1); % Added

% for k = 1:num_density_matrices
%     % ... code for generating rho and correlation_matrix remains unchanged ...
% 
%     % Calculate singular value decomposition (SVD)
%     [U, S, V] = svd(correlation_matrix);
% 
%     % Extract singular values
%     singular_values = diag(S);
% 
%     % Sort singular values in descending order
%     [singular_values_sorted, idx] = sort(singular_values, 'descend');
% 
% %     % Calculate the square root of the sum of squares of the first 100 singular values
% %     iato = norm(singular_values_sorted(1:100));
% %     iato_values(k) = iato; % Save result
% 
%     % ... your original s1, s2, ..., s100 assignments can be omitted
% end

    % Matrix M10
     M10 =[s1 s2 s4 s6 s8 s10 s12 s14 s16 s18; 
      s3 s20 s21 s23 s25 s27 s29 s31 s33 s35; 
      s5 s22 s37 s38 s40 s42 s44 s46 s48 s50; 
      s7 s24 s39 s52 s53 s55 s57 s59 s61 s63; 
      s9 s26 s41 s54 s65 s66 s68 s70 s72 s74;
      s11 s28 s43 s56 s67 s76 s77 s79 s81 s83;
      s13 s30 s45 s58 s69 s78 s85 s86 s88 s90;
      s15 s32 s47 s60 s71 s80 s87 s92 s93 s95;
      s17 s34 s49 s62 s73 s82 s89 s94 s97 s98;
      s19 s36 s51 s64 s75 s84 s91 s96 s99 s100];

    % Calculate eigenvalues
    eigenvalues = eig(M10);

    % Find the maximum eigenvalue
    max_eigenvalue = max(eigenvalues);
    lambda_11 = max_eigenvalue;
    theta_11 = s1+s2+s4+s6+s8+s10+s12+s14+s16+s18-s100*(1-(2*(s19+s36+s51+s64+s75+s84+s91+s96+s99)/(s1+s2+s4+s6+s8+s10+s12+s14+s16+s18-2*s100+sqrt((s1+s2+s4+s6+s8+s10+s12+s14+s16+s18)^2-4*s100*((s1+s2+s4+s6+s8+s10+s12+s14+s16+s18)-(s19+s36+s51+s64+s75+s84+s91+s96+s99+s100))))));
    %zeta_11 = s1+s2+s4+s6+s8+s10+s12+s14+s16+s18-s100*(1-sqrt((s19+s36+s51+s64+s75+s84+s91+s96+s99)/(s1+s2+s4+s6+s8+s10+s12+s14+s16+s18-s100)));
    eta_11 = s1+s2+s4+s6+s8+s10+s12+s14+s16+s18-s100*(1-sqrt((s19+s36+s51+s64+s75+s84+s91+s96+s99+s100)/(s1+s2+s4+s6+s8+s10+s12+s14+s16+s18)));
    P_11 = s1+s2+s4+s6+s8+s10+s12+s14+s16+s18;

    % Store function values
    lambda_11_values(k) = lambda_11;
    theta_11_values(k) = theta_11;
    %zeta_11_values(k) = zeta_11;
    eta_11_values(k) = eta_11;
    P_11_values(k) = P_11;
end

% % Plot function values 1
% figure;
% plot(1:num_density_matrices, lambda_11_values, 'o-', 'DisplayName', '\lambda_11');
% hold on;
% plot(1:num_density_matrices, theta_11_values, 'x-', 'DisplayName', '\theta_11');
% %plot(1:num_density_matrices, zeta_6_values, 's-', 'DisplayName', '\zeta_6');
% plot(1:num_density_matrices, eta_11_values, 'd-', 'DisplayName', '\eta_11');
% plot(1:num_density_matrices, P_11_values, '*-', 'DisplayName', 'P_11');
% xlabel('Density Matrix Index');
% ylabel('Function Value');
% title('Function Values for Random Density Matrices');
% legend show;
% grid on;

% % Plot function values 2
% figure;
% plot(1:num_density_matrices, lambda_11_values, 'o-', 'LineWidth', 1.5, 'MarkerSize', 11, 'DisplayName', '\lambda_11', 'Color', [0 0.4470 0.7410]);
% hold on;
% plot(1:num_density_matrices, theta_11_values, 'x--', 'LineWidth', 1.5, 'MarkerSize', 11, 'DisplayName', '\theta_11', 'Color', [0.8500 0.3250 0.0980]);
% plot(1:num_density_matrices, zeta_11_values, 's-.', 'LineWidth', 1.5, 'MarkerSize', 11, 'DisplayName', '\zeta_11', 'Color', [0.9290 0.6940 0.1250]);
% plot(1:num_density_matrices, eta_11_values, 'd:', 'LineWidth', 1.5, 'MarkerSize', 11, 'DisplayName', '\eta_11', 'Color', [0.4940 0.1840 0.5560]);
% plot(1:num_density_matrices, P_11_values, '*-', 'LineWidth', 1.5, 'MarkerSize', 11, 'DisplayName', 'P_11', 'Color', [0.4660 0.6740 0.1880]);
% xlabel('Density Matrix Index');
% ylabel('Function Value');
% title('Function Values for Random Density Matrices');
% legend show;
% grid on;

% Plot function values 3
figure;
plot(1:num_density_matrices, lambda_11_values, 'o-', 'LineWidth', 1.5, 'MarkerSize', 11, 'DisplayName', '\lambda_{11}', 'Color', [0 0.4470 0.7410]);
hold on;
plot(1:num_density_matrices, theta_11_values, 'x--', 'LineWidth', 1.5, 'MarkerSize', 11, 'DisplayName', '\theta_{11}', 'Color', [0.8500 0.3250 0.0980]);
%plot(1:num_density_matrices, zeta_11_values, 's-.', 'LineWidth', 1.5, 'MarkerSize', 11, 'DisplayName', '\zeta_{11}', 'Color', [0.9290 0.6940 0.1250]);
plot(1:num_density_matrices, eta_11_values, 'd:', 'LineWidth', 1.5, 'MarkerSize', 11, 'DisplayName', '\eta_{11}', 'Color', [0.4940 0.1840 0.5560]);
plot(1:num_density_matrices, P_11_values, '*-', 'LineWidth', 1.5, 'MarkerSize', 11, 'DisplayName', 'P_{11}', 'Color', [0.4660 0.6740 0.1880]);
% Plot iota_5
%plot(1:num_density_matrices, iota_values, 'v-.', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', '\iota_5', 'Color', [0.3010 0.7450 0.9330]);
xlabel('Density Matrix Index');
ylabel('Coefficient and Bounds Values');
title('The optimal witnesses coefficients \lambda_{11} and upper bounds');
legend show;
grid on;

% % Set y-axis range and ticks
% ylim([min([theta_11_values; zeta_11_values; eta_11_values; P_11_values]) - 0.1, max([theta_11_values; zeta_11_values; eta_11_values; P_11_values]) + 0.1]);
% yticks(linspace(min([theta_11_values; zeta_11_values; eta_11_values; P_11_values]) - 0.1, max([theta_11_values; zeta_11_values; eta_11_values; P_11_values]) + 0.1, 10));
% 
% 
% % Calculate the minimum and maximum of all function values
% min_value = min([theta_11_values; zeta_11_values; eta_11_values; P_11_values]);
% max_value = max([theta_11_values; zeta_11_values; eta_11_values; P_11_values]);
% 
% % Set a smaller y-axis range and denser ticks
% ylim([min_value - 0.05, max_value + 0.05]);
% yticks(linspace(min_value - 0.05, max_value + 0.05, 20));

% figure;
% hold on
% plot(1:num_density_matrices, lambda_11_values, 'b-o', 'LineWidth', 1.5, ...
%      'MarkerSize', 6, 'DisplayName', '\lambda_{11}');
% plot(1:num_density_matrices, theta_11_values, 'r--^', 'LineWidth', 1.5, ...
%      'MarkerSize', 6, 'DisplayName', '\theta_{11}');
% plot(1:num_density_matrices, eta_11_values,   'c:*',  'LineWidth', 1.5, ...
%      'MarkerSize', 6, 'DisplayName', '\eta_{11}');
% plot(1:num_density_matrices, P_11_values,     'm-d',  'LineWidth', 1.5, ...
%      'MarkerSize', 6, 'DisplayName', 'P_{11}');
% hold off
% 
% %------------------------- Other original property settings -----------------------------
% title('Optimal witness coefficient \lambda_{11} and related bounds', ...
%       'FontSize',14,'FontWeight','bold');
% xlabel('Density Matrix Index (k)','FontSize',12);
% ylabel('Coefficient and Bounds Values','FontSize',12);
% set(gca,'YScale',y_axis_scale,'FontSize',10);
% 
% % ...... your original y-limit check code goes here ......
% 
% %--------------- Main legend + extra legend in top right ----------------
% hMainLgd  = legend('show','Location','best','FontSize',10);
% hExtraLgd = copyobj(hMainLgd,gcf);
% set(hExtraLgd,'Units','normalized','Position',[0.73 0.77 0.20 0.18]);
% 
% grid on;  set(gca,'GridLineStyle','--');
% 
% %--------------- Compress vertical spacing (choose one) ----------------
% % Method 1: Compress the axes height
% axPos = get(gca,'Position');  axPos(4) = axPos(4)*0.60;  set(gca,'Position',axPos);

% % Generate a random 100x100 matrix
% correlation_matrix = randn(100);
% 
% % Calculate singular value decomposition (SVD)
% [U, S, V] = svd(correlation_matrix);
% 
% % Extract singular values
% singular_values = diag(S);
% 
% % Sort singular values in descending order
% [singular_values_sorted, idx] = sort(singular_values, 'descend');
% 
% % Assign the first 100 singular values to s1, s2, ..., s100
% for i = 1:100
%     eval(['s' num2str(i) ' = singular_values_sorted(i);']);
% end
