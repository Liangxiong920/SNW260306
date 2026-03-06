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
% matrix. The optimal witness coefficient \lambda_{16} and upper bounds 
% of the mixed state are derived and analyzed.
% =========================================================================



% Set random seed for reproducibility
rng(1);

% Define system dimensions
m = 15;
n = 15;

% Define the number of density matrices to generate
num_density_matrices = 50;

% Initialize arrays to store function values for each density matrix
lambda_16_values = zeros(num_density_matrices, 1);
theta_16_values = zeros(num_density_matrices, 1);
%zeta_16_values = zeros(num_density_matrices, 1);
eta_16_values = zeros(num_density_matrices, 1);
P_16_values = zeros(num_density_matrices, 1);

for k = 1:num_density_matrices
    % Initialize the density matrix
    rho = zeros(m*n, m*n);

    % Define the number of pure states to generate
    num_pure_states = 10000;

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
    rho_vector = transpose(rho(:));

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

    % Assign the first 225 singular values to s1, s2, ..., s225
    for i = 1:225
        eval(['s' num2str(i) ' = singular_values_sorted(i);']);
    end

    % Define M15 matrix
    M15 = [ ...
        s1   s2   s4   s6   s8   s10  s12  s14  s16  s18  s20  s22  s24  s26  s28 ; ...
        s3   s30  s31  s33  s35  s37  s39  s41  s43  s45  s47  s49  s51  s53  s55 ; ...
        s5   s32  s57  s58  s60  s62  s64  s66  s68  s70  s72  s74  s76  s78  s80 ; ...
        s7   s34  s59  s82  s83  s85  s87  s89  s91  s93  s95  s97  s99  s101 s103; ...
        s9   s36  s61  s84  s105 s106 s108 s110 s112 s114 s116 s118 s120 s122 s124; ...
        s11  s38  s63  s86  s107 s126 s127 s129 s131 s133 s135 s137 s139 s141 s143; ...
        s13  s40  s65  s88  s109 s128 s145 s146 s148 s150 s152 s154 s156 s158 s160; ...
        s15  s42  s67  s90  s111 s130 s147 s162 s163 s165 s167 s169 s171 s173 s175; ...
        s17  s44  s69  s92  s113 s132 s149 s164 s177 s178 s180 s182 s184 s186 s188; ...
        s19  s46  s71  s94  s115 s134 s151 s166 s179 s190 s191 s193 s195 s197 s199; ...
        s21  s48  s73  s96  s117 s136 s153 s168 s181 s192 s201 s202 s204 s206 s208; ...
        s23  s50  s75  s98  s119 s138 s155 s170 s183 s194 s203 s210 s211 s213 s215; ...
        s25  s52  s77  s100 s121 s140 s157 s172 s185 s196 s205 s212 s217 s218 s220; ...
        s27  s54  s79  s102 s123 s142 s159 s174 s187 s198 s207 s214 s219 s222 s223; ...
        s29  s56  s81  s104 s125 s144 s161 s176 s189 s200 s209 s216 s221 s224 s225 ...
    ];

    % Calculate eigenvalues
    eigenvalues = eig(M15);
    P_16 = s1 + s2 + s4 + s6 + s8 + s10 + s12 + s14 + s16 + s18 + s20 + s22 + s24 + s26 + s28;
    p_16 = s29 + s56 + s81 + s104 + s125 + s144 + s161 + s176 + s189 + s200 + s209 + s216 + s221 + s224 + s225;


    % Find the maximum eigenvalue
    max_eigenvalue = max(eigenvalues);
    lambda_16 = max_eigenvalue;

    % Calculation expressions for theta_16, zeta_16, eta_16, and P_16
    theta_16 = P_16 - s225 * (1 - (2 * (p_16-s225) / (P_16 - 2 * s225 + sqrt((P_16)^2 - 4 * s225 * ((P_16) - (p_16))))));

    %zeta_16 = P_16 - s225 * (1 - sqrt((p_16-s225) / (P_16 - s225)));

    eta_16 = P_16 - s225 * (1 - sqrt((p_16) / (P_16)));

    P_16 = s1 + s2 + s4 + s6 + s8 + s10 + s12 + s14 + s16 + s18 + s20 + s22 + s24 + s26 + s28;

    % Store function values
    lambda_16_values(k) = lambda_16;
    theta_16_values(k) = theta_16;
    %zeta_16_values(k) = zeta_16;
    eta_16_values(k) = eta_16;
    P_16_values(k) = P_16;
end

% % Calculate the minimum and maximum of all function values to set the y-axis range
% all_values = [lambda_16_values; theta_16_values;  eta_16_values; P_16_values];
% y_min = min(all_values);
% y_max = max(all_values);
% 
% % If the value range is large, consider using a logarithmic scale
% if y_max / y_min > 100
%     y_axis_scale = 'log';
% else
%     y_axis_scale = 'linear';
% end
% 
% % Plot image - Compare the five functions in a single figure
% figure;
% hold on;
% plot(1:num_density_matrices, lambda_16_values, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', '\lambda_{16}');
% plot(1:num_density_matrices, theta_16_values, 'r--^', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', '\theta_{16}');
% %plot(1:num_density_matrices, zeta_16_values, 'g-.s', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', '\zeta_{16}');
% plot(1:num_density_matrices, eta_16_values, 'c:*', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', '\eta_{16}');
% plot(1:num_density_matrices, P_16_values, 'm-d', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'P_{16}');
% hold off;
% 
% % Set figure properties
% title('Optimal witness coefficient \lambda_{16} and related bounds', 'FontSize', 14, 'FontWeight', 'bold');
% xlabel('Density Matrix Index (k)', 'FontSize', 12);
% ylabel('Coefficient and Bounds Values', 'FontSize', 12);
% set(gca, 'YScale', y_axis_scale, 'FontSize', 10);
% if strcmp(y_axis_scale, 'linear')
%     %ylim([y_min - 0.05*abs(y_min), y_max + 0.05*abs(y_max)]); % Add 5% margin
%     y_min = min(dataCell{k});
%     y_max = max(dataCell{k});
% if y_max == y_min
%     ylim([y_min-1e-6, y_max+1e-6]);
%     else
%     ylim([y_min - 0.05*abs(y_min), y_max + 0.05*abs(y_max)]);
%     end
% end
% legend('show', 'Location', 'best', 'FontSize', 10);
% grid on;
% set(gca, 'GridLineStyle', '--');


figure;
hold on
plot(1:num_density_matrices, lambda_16_values, 'b-o', 'LineWidth', 1.5, ...
     'MarkerSize', 6, 'DisplayName', '\lambda_{16}');
plot(1:num_density_matrices, theta_16_values, 'r--^', 'LineWidth', 1.5, ...
     'MarkerSize', 6, 'DisplayName', '\theta_{16}');
plot(1:num_density_matrices, eta_16_values,   'c:*',  'LineWidth', 1.5, ...
     'MarkerSize', 6, 'DisplayName', '\eta_{16}');
plot(1:num_density_matrices, P_16_values,     'm-d',  'LineWidth', 1.5, ...
     'MarkerSize', 6, 'DisplayName', 'P_{16}');
hold off

%------------------------- Other original property settings -----------------------------
title('Optimal witness coefficient \lambda_{16} and related bounds', ...
      'FontSize',14,'FontWeight','bold');
xlabel('Density Matrix Index (k)','FontSize',12);
ylabel('Coefficient and Bounds Values','FontSize',12);
set(gca,'YScale',y_axis_scale,'FontSize',10);

% ...... your original y-limit check code goes here ......

%--------------- Main legend + extra legend in top right ----------------
hMainLgd  = legend('show','Location','best','FontSize',10);
hExtraLgd = copyobj(hMainLgd,gcf);
set(hExtraLgd,'Units','normalized','Position',[0.73 0.77 0.20 0.18]);

grid on;  set(gca,'GridLineStyle','--');

%--------------- Compress vertical spacing (choose one) ----------------
% Method 1: Compress the axes height
axPos = get(gca,'Position');  axPos(4) = axPos(4)*0.60;  set(gca,'Position',axPos);

% Method 2 (use when Method 1 is commented out): Control vertical aspect ratio
% set(gca,'PlotBoxAspectRatio',[1 0.30 1]);