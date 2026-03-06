% =========================================================================
% NUMERICAL INVESTIGATION OF RANDOM MIXED STATES AND WITNESS BOUNDS
%
% This script numerically investigates different mixed states by generating 
% 50 independent samples using MATLAB. Each random density matrix is 
% constructed as a finite convex combination of independently sampled 
% Haar-distributed pure states. Specifically, for each sample, 50 independent 
% Haar-random pure states |\psi_j> are generated via QR decomposition of 
% complex Ginibre matrices. Random positive weights are assigned and the 
% resulting operator is normalized to unit trace, ensuring a valid density 
% matrix. The optimal witness coefficient \lambda_{21} and upper bounds 
% of the mixed state are derived and analyzed.
% =========================================================================

% Set random seed for reproducibility
rng(1);

% Define system dimensions
m = 20;
n = 20;

% Define the number of density matrices to generate
num_density_matrices = 50;

% Initialize arrays to store function values for each density matrix
lambda_20_values = zeros(num_density_matrices, 1);
theta_20_values = zeros(num_density_matrices, 1);
%zeta_20_values = zeros(num_density_matrices, 1);
eta_20_values = zeros(num_density_matrices, 1);
P_20_values = zeros(num_density_matrices, 1);

for k = 1:num_density_matrices
    % Initialize the density matrix
    rho = zeros(m*n, m*n);

    % Define the number of pure states to generate
    num_pure_states = 50;

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

    % Assign the first 400 singular values to s1, s2, ..., s400
    for i = 1:400
        eval(['s' num2str(i) ' = singular_values_sorted(i);']);
    end

  
    M20 = [ ...
    s1   s2   s4   s6   s8   s10  s12  s14  s16  s18  s20  s22  s24  s26  s28  s30  s32  s34  s36  s38 ; ...
    s3   s40  s41  s43  s45  s47  s49  s51  s53  s55  s57  s59  s61  s63  s65  s67  s69  s71  s73  s75 ; ...
    s5   s42  s77  s78  s80  s82  s84  s86  s88  s90  s92  s94  s96  s98  s100 s102 s104 s106 s108 s110; ...
    s7   s44  s79  s112 s113 s115 s117 s119 s121 s123 s125 s127 s129 s131 s133 s135 s137 s139 s141 s143; ...
    s9   s46  s81  s114 s145 s146 s148 s150 s152 s154 s156 s158 s160 s162 s164 s166 s168 s170 s172 s174; ...
    s11  s48  s83  s116 s147 s176 s177 s179 s181 s183 s185 s187 s189 s191 s193 s195 s197 s199 s201 s203; ...
    s13  s50  s85  s118 s149 s178 s205 s206 s208 s210 s212 s214 s216 s218 s220 s222 s224 s226 s228 s230; ...
    s15  s52  s87  s120 s151 s180 s207 s232 s233 s235 s237 s239 s241 s243 s245 s247 s249 s251 s253 s255; ...
    s17  s54  s89  s122 s153 s182 s209 s234 s257 s258 s260 s262 s264 s266 s268 s270 s272 s274 s276 s278; ...
    s19  s56  s91  s124 s155 s184 s211 s236 s259 s280 s281 s283 s285 s287 s289 s291 s293 s295 s297 s299; ...
    s21  s58  s93  s126 s157 s186 s213 s238 s261 s282 s301 s302 s304 s306 s308 s310 s312 s314 s316 s318; ...
    s23  s60  s95  s128 s159 s188 s215 s240 s263 s284 s303 s320 s321 s323 s325 s327 s329 s331 s333 s335; ...
    s25  s62  s97  s130 s161 s190 s217 s242 s265 s286 s305 s322 s337 s338 s340 s342 s344 s346 s348 s350; ...
    s27  s64  s99  s132 s163 s192 s219 s244 s267 s288 s307 s324 s339 s352 s353 s355 s357 s359 s361 s363; ...
    s29  s66  s101 s134 s165 s194 s221 s246 s269 s290 s309 s326 s341 s354 s365 s366 s368 s370 s372 s374; ...
    s31  s68  s103 s136 s167 s196 s223 s248 s271 s292 s311 s328 s343 s356 s367 s376 s377 s379 s381 s383; ...
    s33  s70  s105 s138 s169 s198 s225 s250 s273 s294 s313 s330 s345 s358 s369 s378 s385 s386 s388 s390; ...
    s35  s72  s107 s140 s171 s200 s227 s252 s275 s296 s315 s332 s347 s360 s371 s380 s387 s392 s393 s395; ...
    s37  s74  s109 s142 s173 s202 s229 s254 s277 s298 s317 s334 s349 s362 s373 s382 s389 s394 s397 s398; ...
    s39  s76  s111 s144 s175 s204 s231 s256 s279 s300 s319 s336 s351 s364 s375 s384 s391 s396 s399 s400 ...
];

    % Calculate eigenvalues and related sums
    eigenvalues = eig(M20);

    P_21=s1 + s2 + s4 + s6 + s8 + s10 + s12 + s14 + s16 + s18 + s20 + s22 + s24 + s26 + s28+s30+s32+s34+s36+s38;

    p_21=s39+s76+s111+s144+s175+s204+s231+s256+s279+s300+s319+s336+s351+s364+s375+s384+s391+s396+s399+s400;

    % Find the maximum eigenvalue
    max_eigenvalue = max(eigenvalues);
    lambda_21 = max_eigenvalue;

    % Calculation expressions for theta_21, zeta_21, eta_21, and P_21
    theta_21 = P_21 - s400 * (1 - (2 * (p_21-s400) / (P_21 - 2 * s400 + sqrt((P_21)^2 - 4 * s400 * ((P_21) - (p_21))))));

    %zeta_21 = P_21 - s400 * (1 - sqrt((p_21-s400) / (P_21-s400)));

    eta_21 = P_21 - s400 * (1 - sqrt((p_21) / (P_21)));

    P_21 = s1 + s2 + s4 + s6 + s8 + s10 + s12 + s14 + s16 + s18 + s20 + s22 + s24 + s26 + s28+s30+s32+s34+s36+s38;

    % Store function values
    lambda_21_values(k) = lambda_21;
    theta_21_values(k) = theta_21;
    %zeta_21_values(k) = zeta_21;
    eta_21_values(k) = eta_21;
    P_21_values(k) = P_21;
end



% % Plot image - Compare the functions in a single figure
% figure;
% hold on;
% plot(1:num_density_matrices, lambda_21_values, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', '\lambda_{21}');
% plot(1:num_density_matrices, theta_21_values, 'r--^', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', '\theta_{21}');
% plot(1:num_density_matrices, eta_21_values, 'c:*', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', '\eta_{21}');
% plot(1:num_density_matrices, P_21_values, 'm-d', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'P_{21}');
% hold off;
% 
% % Set figure properties
% title('Optimal witness coefficient \lambda_{21} and related bounds', 'FontSize', 14, 'FontWeight', 'bold');
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
plot(1:num_density_matrices, lambda_21_values, 'b-o', 'LineWidth', 1.5, ...
     'MarkerSize', 6, 'DisplayName', '\lambda_{21}');
plot(1:num_density_matrices, theta_21_values, 'r--^', 'LineWidth', 1.5, ...
     'MarkerSize', 6, 'DisplayName', '\theta_{21}');
plot(1:num_density_matrices, eta_21_values,   'c:*',  'LineWidth', 1.5, ...
     'MarkerSize', 6, 'DisplayName', '\eta_{21}');
plot(1:num_density_matrices, P_21_values,     'm-d',  'LineWidth', 1.5, ...
     'MarkerSize', 6, 'DisplayName', 'P_{21}');
hold off

%------------------------- Other original property settings -----------------------------
title('Optimal witness coefficient \lambda_{21} and related bounds', ...
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