function [A_matrix, B_matrix, C_matrix, D_matrix, Voltage] = Sparse_state_space_solution(Data)
%{
Optimized function for solving the state-space model of the cascaded Pi circuit using sparse matrices
and the trapezoidal rule with LU decomposition.

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
V_modal = Data.V       % Voltage in the modal domain, ab0

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

	%                         Input voltage source                            %
	% ─────────────────────────────────────────────────────────────────────────
	u = V_modal;   % Voltage source
	% ─────────────────────────────────────────────────────────────────────────


	%                         State-space submatrices                          %
	% ─────────────────────────────────────────────────────────────────────────

	%                            Identity matrices                            %
	I   = speye(M);


	%                            Submatrices of A                             %
	% A_11
	A_11 = sparse(-R_in/L_in);

	% A_12
	A_12 = sparse(1,M);

	% A_13
	A_13 = sparse(1,M*m);

	% A_14
	A_14 = sparse(1,M+1);
	A_14(1,1) = -1/L_in;

	% A_21
	A_21 = sparse(M,1);

	% A_31
	A_31 = sparse(M*m,1);

	% A_41
	A_41 = sparse(M+1,1);
	A_41(1,1) = 1/C_in;

	% A_34
	A_34 = sparse(M*m,M+1);

	% A_43
	A_43 = sparse(M+1,M*m);

	% A_22
	R_sum = -(sum(Rk)+Ro);
	A_22 = (R_sum/Lo)*I;

	% A_32
	A_32 = sparse(M*m,M);
	for k = 1:m
		rows = (k-1)*M + (1:M);
		A_32(rows,:) = (Rk(k)/Lk(k))*I;
	end

	% A_23
	A_23 = sparse(M,M*m);
	for k = 1:m
		cols = (k-1)*M + (1:M);
		A_23(:,cols) = (Rk(k)/Lo)*I;
	end

	% A_44
	if M == 1
		A_44 = sparse(2,2);
		A_44(1,1) = -G_in/C_in;
		A_44(2,2) = -(R_L*G_fin + 1)/(R_L*C_fin);
	else
		A_44 = sparse(M+1,M+1);
		A_44(1,1) = -G_in/C_in;
		for k = 2:M
			A_44(k,k) = -G/Ck;
		end
		A_44(M+1,M+1) = -(R_L*G_fin + 1)/(R_L*C_fin);
	end

	% A_33
	A_33 = sparse(M*m,M*m);
	for k = 1:m
		idx = (k-1)*M + (1:M);
		A_33(idx,idx) = (-Rk(k)/Lk(k))*I;
	end

	% A_42
	A_42 = sparse(M+1,M);
	for k = 1:M
		if k == 1
			A_42(k,k)   = -1/C_in;
			A_42(k+1,k) =  1/Ck;
		elseif k == M
			A_42(k,k)   = -1/Ck;
			A_42(k+1,k) =  1/C_fin;
		else
			A_42(k,k)   = -1/Ck;
			A_42(k+1,k) =  1/Ck;
		end
	end

	% A_24
	A_24 = sparse(M,M+1);
    for k = 1:M
		A_24(k,k)   =  1/Lo;
		A_24(k,k+1) = -1/Lo;
    end
    % ─────────────────────────────────────────────────────────────────────────


	%                           State-space matrices                         %
	% ─────────────────────────────────────────────────────────────────────────
	% Matrix A
	A = [A_11  A_12  A_13  A_14;
		 A_21  A_22  A_23  A_24;
		 A_31  A_32  A_33  A_34;
		 A_41  A_42  A_43  A_44];
	A_size = size(A,1);
	A_matrix.A = A;

	% Matrix B
	B = sparse(A_size,1); B(1,1) = 1/L_in; B_matrix.B = B;

	% Matrix C
	C_mat = sparse(1,A_size); C_mat(end) = 1; C_matrix.C = C_mat;

	% Matrix D
	D = 0; D_matrix.D = D;
	% ─────────────────────────────────────────────────────────────────────────


	%                             Trapezoidal-rule                            %
	% ─────────────────────────────────────────────────────────────────────────
	Mmat    = speye(A_size) - (dt/2)*A;
	Nmat = speye(A_size) + (dt/2)*A;

	% LU factorization
	[L,U,P,Q] = lu(Mmat);

	% Variable initialization
	x_trap = zeros(A_size,N); % Initial conditions
	y = zeros(1,N);           % Outputs

	% Main loop
	for k = 1:N-1
		rhs = Nmat*x_trap(:,k) + (dt/2)*B*(u(:,k)+u(:,k+1));
		x_trap(:,k+1) = Q*(U\(L\(P*rhs)));
		y(:,k) = C_mat*x_trap(:,k);
	end
	y(:,N) = C_mat*x_trap(:,N);
	Voltage = y;
	% ─────────────────────────────────────────────────────────────────────────


	%                      End of the simulation timer                       %
	% ─────────────────────────────────────────────────────────────────────────
	simulation_time = toc;
	fprintf('Sparse state-space simulation completed in %.4f seconds\n', simulation_time);
	% ─────────────────────────────────────────────────────────────────────────


	%                                  Plots                                 %
	% ─────────────────────────────────────────────────────────────────────────  
	% Sparsity plot
	if M <= 5 && m <= 5
		figure;
		spy(A,'k.',40)
		xlabel('Columns');
		ylabel('Rows');
		set(gcf,'pos',[200 200 800 800]);
	end
	% ───────────────────────────────────────────────────────────────────────── 
end
