function [A_matrix, B_matrix, C_matrix, D_matrix, Voltage] = State_space_solution(Data)
%{
Function that solves the state-space model of the cascaded Pi circuit
for the single-phase transmission line

Inputs:
Data - Structure with the following parameters:
N       = Data.N;      % Number of samples
dt      = Data.dt;     % Time step
R_in    = Data.R_in;   % Input resistance
Ro      = Data.Ro;     % Series resistance
R_k     = Data.Rk;     % Parallel resistance
R_L     = Data.R_L;    % Load resistance
L_in    = Data.L_in;   % Input inductance
Lo      = Data.Lo;     % Series inductance
L_k     = Data.Lk;     % Parallel inductance
Ck      = Data.Ck;     % Capacitance
G       = Data.G;      % Conductance
M       = Data.M;      % Number of cascaded Pi circuits
V_modal = Data.V       % Voltage in the modal domanin, ab0

Outputs:
A    - State-space matrices
V    - Voltage in R_L
%}


	%                     Start of the simulation timer                      %
	% ─────────────────────────────────────────────────────────────────────────
	tic;
	% ─────────────────────────────────────────────────────────────────────────


	%                         Parameter extraction                           %
	% ─────────────────────────────────────────────────────────────────────────
	N        = Data.N;    % Number of samples
	dt       = Data.dt;   % Time step
	R_in     = Data.R_in; % Input resistance
	Ro       = Data.Ro;   % Series resistance
	Rk       = Data.Rk;   % Parallel resistance
	R_L      = Data.R_L;  % Load resistance
	L_in     = Data.L_in; % Input inductance
	Lo       = Data.Lo;   % Series inductance
	Lk       = Data.Lk;   % Parallel inductance
	Ck       = Data.Ck;   % Capacitance
	G        = Data.G;    % Conductance
	M        = Data.M;    % Number of cascaded Pi circuits
	V_modal  = Data.V;    % Voltage in the modal domain, ab0
	% ─────────────────────────────────────────────────────────────────────────


	%                               Auxiliary data                            %
	% ─────────────────────────────────────────────────────────────────────────
	% Number of fitting poles (m)
	m = length(Rk);
	% Initial and final capacitances and conductances
	C_in = Ck/2;   C_fin = C_in;   G_in = G/2;   G_fin = G_in;
	% ─────────────────────────────────────────────────────────────────────────


	%                         Input voltage source                          %
	% ─────────────────────────────────────────────────────────────────────────
	u = V_modal;   % Voltage source
	% ─────────────────────────────────────────────────────────────────────────


	%                         State-space submatrices                        %
	% ─────────────────────────────────────────────────────────────────────────
	% Auxiliary matrices
	I   = eye(M);
	I_1 = eye(M-1);
	
	% Submatrices of A
	A_11 = -R_in/L_in;
	A_12 = zeros(1, M);
	A_13 = zeros(1, M*m);
	A_14 = [-1/L_in, zeros(1,M)];
	A_21 = zeros(M,1);
	A_31 = zeros(M*m, 1);
	A_41 = [1/C_in; zeros(M,1)];
	A_34 = zeros(M*m, M+1);
	A_43 = zeros(M+1, M*m);
	
	R_sum = -(sum(Rk)+Ro);
	A_22 = (R_sum/Lo) * I;
	
	% A_32
	A_32 = zeros(M*m, M);
	for k = 1:m
		start_row = (k-1)*M + 1;
		end_row = M*k;
		A_aux = I * (Rk(k)/Lk(k));
		A_32(start_row:end_row, :) = A_aux;
	end
	
	% A_23
	A_23 = zeros(M, M*m);
	for k = 1:m
		start_col = (k-1)*M + 1;
		end_col = M*k;
		A_aux = I * (Rk(k)/Lo);
		A_23(:, start_col:end_col) = A_aux;
	end
	
	% A_44
	A_44 = [-G_in/C_in, zeros(1, M-1), 0;
			zeros(M-1, 1), I_1*(-G/Ck), zeros(M-1, 1);
			0, zeros(1, M-1), -(R_L*G_fin + 1)/(R_L * C_fin)];
	if M == 1
		A_44 = [-G_in/C_in, 0;
				0, -(R_L*G_fin + 1)/(R_L * C_fin)];
	end
	
	% A_33
	A_33 = zeros(M*m, M*m);
	for l = 1:m
		start_idx = (l-1)*M + 1;
		end_idx = M*l;
		A_aux = I * (-Rk(l)/Lk(l));
		A_33(start_idx:end_idx, start_idx:end_idx) = A_aux;
	end
	
	% A_42
	A_42 = zeros(M+1, M);
	for k = 1:M
		start_idx = k;
		end_idx = k+1;
		if k == 1
			A_aux = [-1/C_in; 1/Ck];
		elseif k == M
			A_aux = [-1/Ck; 1/C_fin];
		else
			A_aux = [-1/Ck; 1/Ck];
		end
		A_42(start_idx:end_idx, k) = A_aux;
	end
	if M == 1
		A_aux = [-1/C_in; 1/C_fin];
		A_42 = A_aux;
	end
	
	% A_24
	A_24 = zeros(M, M+1);
	for k = 1:M
		start_idx = k;
		end_idx = k+1;
		A_aux = [1/Lo, -1/Lo];
		A_24(k, start_idx:end_idx) = A_aux;
	end
	% ─────────────────────────────────────────────────────────────────────────


	%                           State-space matrices                         %
	% ─────────────────────────────────────────────────────────────────────────
	% Matrix A
	A = [A_11, A_12, A_13, A_14;
		 A_21, A_22, A_23, A_24;
		 A_31, A_32, A_33, A_34;
		 A_41, A_42, A_43, A_44];

	A_size = length(A);

	A_matrix.A = A; % Structure
	
	% Matrix B
	B = [1/L_in;
		 zeros(M,1);
		 zeros(M*m,1);
		 zeros(M+1,1)];
	B_matrix.B = B; % Structure

	% Matrix C
	C_mat = [zeros(1, A_size-1), 1];
	C_matrix.C = C_mat; % Structure

	% Matrix D
	D = 0;
	D_matrix.D = D; % Structure
	% ─────────────────────────────────────────────────────────────────────────


	%                            Trapezoidal-rule                             %
	% ─────────────────────────────────────────────────────────────────────────
	I_mat = eye(size(A));
	Alpha = (I_mat - dt/2 * A) \ (I_mat + dt/2 * A);
	Beta = (I_mat - dt/2 * A) \ (dt/2 * B);

	% Variable initialization
	x_trap = zeros(length(A), N);
	y = zeros(1, N);

	% Iterations
	for k = 1:N-1
		x_trap(:, k+1) = Alpha * x_trap(:, k) + Beta * (u(:,k) + u(:,k+1));
		y(:,k) = C_mat * x_trap(:, k) + D * u(:,k);
	end
	y(:,N) = C_mat * x_trap(:, N) + D * u(:,N);

	Voltage = y;
	% ─────────────────────────────────────────────────────────────────────────


	%                      End of the simulation timer                       %
	% ─────────────────────────────────────────────────────────────────────────
	simulation_time = toc;
	fprintf('State-space simulation completed in:       %.4f seconds\n', simulation_time);
	% ─────────────────────────────────────────────────────────────────────────


	%                                  Plots                                 %
	% ─────────────────────────────────────────────────────────────────────────  
	% Sparsity of A
	if M <= 5 && m <= 5
		figure;
		spy(A, 'k.', 40)
		xlabel('Columns'); ylabel('Rows');
		set(gcf,'pos',[200 200 800 800]);
	end
	% ─────────────────────────────────────────────────────────────────────────

end