function [Gm_matrix,V_Rl,c_ind,c_cap] = EMT_solution(Data)

% This function solves a single-phase transmission line in the time domain
% represented by cascaded Pi circuits using an EMT-type formulation.

% Inputs:
%     Data - Structure with the following parameters:
%         N             = Data.N;             % Number of samples
%         dt            = Data.dt;            % Time step
%         R_in          = Data.R_in;          % Input resistance
%         Ro            = Data.Ro;            % Series resistance R0
%         Rk            = Data.Rk;            % Rk resistances from the rational fitting
%         R_L           = Data.R_L;           % Load resistance
%         L_in          = Data.L_in;          % Input inductance
%         Lo            = Data.Lo;            % Series inductance L0
%         Lk            = Data.Lk;            % Lk inductances from the rational fitting
%         Ck            = Data.Ck;            % Capacitance of each Pi section
%         n             = Data.n;             % Number of cascaded Pi circuits
%         V             = Data.V;             % Input voltage in the modal domain
%         solution_type = Data.solution_type; % 1.- Thomas    2.- Inverse of Gm
%
% Outputs:
%     Gm_matrix - Structure containing the conductance matrix Gm
%     V_Rl      - Voltage across the load resistance R_L
%     c_ind     - Total number of inductors in the equivalent circuit
%     c_cap     - Total number of capacitors in the equivalent circuit


	%                     Start of the simulation timer                       %
	% ─────────────────────────────────────────────────────────────────────────
	tic;
	% ─────────────────────────────────────────────────────────────────────────


	%                         Parameter extraction                            %
	% ─────────────────────────────────────────────────────────────────────────
	N       = Data.N;     % Number of samples
	dt      = Data.dt;    % Time step
	R_in    = Data.R_in;  % Input resistance
	Ro      = Data.Ro;    % Series resistance
	Rk      = Data.Rk;    % Parallel resistance
	R_L     = Data.R_L;   % Load resistance
	L_in    = Data.L_in;  % Input inductance
	Lo      = Data.Lo;    % Series inductance
	Lk      = Data.Lk;    % Parallel inductance
	Ck      = Data.Ck;    % Capacitance
	%G      = Data.G;     % Conductance
	n_pi    = Data.n;     % Number of cascaded Pi circuits
	V_modal = Data.V;     % Voltage in the ab0 modal domain
	solution_type = Data.solution_type; % 1.- Thomas 2.- Inverse of Gm
	% ─────────────────────────────────────────────────────────────────────────


	%                         Input voltage source                            %
	% ─────────────────────────────────────────────────────────────────────────
	V_source = V_modal;
	% ─────────────────────────────────────────────────────────────────────────


	%                       Conductance-matrix parameters                     %
	% ─────────────────────────────────────────────────────────────────────────
	Np = length(Rk);                   % Number of fitting poles
	Nk = Np - 1;                       % Number of nodes in the RL_k branches
	num_nodes = 2 + n_pi * (3 + Nk);   % Total number of nodes

	% R, L, and C parameters of the equivalent circuit
	Rt = [Ro; Rk];       % Fitting resistances
	Lt = [Lo; Lk];       % Fitting inductances
	Cin = Ck/2;          % Capacitance inicial
	Cfin = Cin;          % Capacitance final

	% Equivalent conductances of R, L, and C
	Gr_in = 1/R_in;      % Input-resistance conductance
	GL_in = dt/(2*L_in); % Equivalent conductance of L_in
	Grl = 1/R_L;         % Resistive-load conductance
	Grt = 1./Rt;         % Conductances of the fitting resistances
	GLt = dt./(2*Lt);    % Equivalent conductances of the inductors
	GC = 2*Ck/dt;        % Equivalent conductance of C
	GC_in = 2*Cin/dt;    % Equivalent conductance of C_in
	GC_fin = 2*Cfin/dt;  % Equivalent conductance of C_fin
	% ─────────────────────────────────────────────────────────────────────────


	%                      Conductance-matrix assembly                       %
	% ─────────────────────────────────────────────────────────────────────────
	% Equivalent conductances G0, G1, ..., GN
	Gk = zeros(1,Np+1); Gk(1,1) = 1/Ro + dt/(2*Lo);
	for i = 1:Np; Gk(i+1) = 1./Rk(i,1) + dt./(2*Lk(i,1)); end

	% Consecutive sums: G1+G2, G2+G3, ..., G(N-1)+GN
	Gk_aux = zeros(1,Np-1);
	for i = 1:Np-1; Gk_aux(1,i) = Gk(1,i+1) + Gk(1,i+2); end

	% Vector W: upper and lower diagonals of Gm
    % w1
	w1 = [GL_in, Grt(1,1), GLt(1,1), Gk(1,2:Np+1)];

	% wm
	wm = [Grt(1,1), GLt(1,1), Gk(1,2:Np+1)];

	% Complete upper and lower diagonal vector
	W = [w1, repmat(wm,1,n_pi-1)]';

	% Vector Q: main diagonal of Gm
	% q1 
	q1 = - [Gr_in + GL_in, Grt(1,1) + GC_in + GL_in, Gk(1,1), GLt(1,1) + Gk(1,2), Gk_aux, Grt(1,1) + Gk(1,Np+1) + GC];

	% qz
	qz = - [Gk(1,1),  GLt(1,1) + Gk(1,2), Gk_aux, Grt(1,1) + Gk(1,Np+1) + GC];

	% qm
	qm = - [Gk(1,1), GLt(1,1) + Gk(1,2), Gk_aux, Gk(1,Np+1) + Grl + GC_fin];

	% Complete main-diagonal vector
	Q = [q1, repmat(qz,1,n_pi-2), qm]';

	% Tridiagonal conductance matrix
	Gm = diag(W,-1) + diag(Q) + diag(W,1);

	Gm_matrix.Gm = Gm;
	% ─────────────────────────────────────────────────────────────────────────


	%                         Precomputation of Gm inverse                    %
	% ─────────────────────────────────────────────────────────────────────────
	if solution_type == 2
		Gm_inv = inv(Gm);
	end
	% ─────────────────────────────────────────────────────────────────────────


	%                           Simulation elements                           %
	% ─────────────────────────────────────────────────────────────────────────
	% Number of inductors and capacitors
	c_ind = n_pi*(Np + 1) + 1; % M(N+1)+1
	c_cap = n_pi + 1;          % M+1

	% Equivalent conductances of all inductors
	GL_tot = [GL_in; repmat(GLt,n_pi,1)];

	% Equivalent conductances of all capacitors
	GC_tot = [GC_in; repmat(GC,n_pi-1,1); GC_fin];

	% General incidence matrix of the series branches
	num_branches = num_nodes - 1;
	rows    = [1:num_branches, 2:num_nodes];
	columns = [1:num_branches, 1:num_branches];
	values  = [ones(1,num_branches), -ones(1,num_branches)];
	B       = sparse(rows,columns,values,num_nodes,num_branches);

	% Inductor incidence matrix
	L_pattern = [true; repmat([false; true(Np+1,1)],n_pi,1)];
	B_L = B(:,L_pattern);

	% Capacitor incidence matrix
	C_nodes = (2 + (0:n_pi)*(Np+2)).';
	B_C = sparse(C_nodes,(1:c_cap).',ones(c_cap,1),num_nodes,c_cap);

	% Inductor and capacitor currents
	L_currents = zeros(N,c_ind);    C_currents = zeros(N,c_cap);
	% ─────────────────────────────────────────────────────────────────────────


	%                         Matrix initialization                           %
	% ─────────────────────────────────────────────────────────────────────────
	% Inductor and capacitor history terms
	L_history = zeros(N,c_ind);     C_history = zeros(N,c_cap);

	% Nodal voltages and history vector
	V = zeros(num_nodes,N);         H = zeros(num_nodes,1);
	% ─────────────────────────────────────────────────────────────────────────


	%                         Nodal-voltage calculation                       %
	% ─────────────────────────────────────────────────────────────────────────
	for k = 2:N

		% Inductor and capacitor history terms
		V_L_prev = B_L.' * V(:,k-1);
		V_C_prev = B_C.' * V(:,k-1);

		% Update of the L and C history terms
		L_history(k,:) = (L_currents(k-1,:).' + GL_tot .* V_L_prev).';
		C_history(k,:) = (-C_currents(k-1,:).' - GC_tot .* V_C_prev).';

		% History vector H
		H_L = B_L * L_history(k,:).';
		H_C = B_C * C_history(k,:).';
		H = H_L + H_C;

		% Contribution of the voltage source at the input node
		H(1,1) = H(1,1) - Gr_in * V_source(1,k);

		% Solution of the nodal voltages
		if solution_type == 1
			V(:,k) = Thomas_Algorithm(Gm,H);
		elseif solution_type == 2
			V(:,k) = Gm_inv*H;
		end

		% Inductor and capacitor currents
		V_L = B_L.' * V(:,k);
		V_C = B_C.' * V(:,k);

		% Update of the L and C currents
		L_currents(k,:) = (L_history(k,:).' + GL_tot .* V_L).';
		C_currents(k,:) = (C_history(k,:).' + GC_tot .* V_C).';

	end

	% Voltage across the load resistance
	V_size = length(V(:,1));
	V_Rl = V(V_size,:);
	% ─────────────────────────────────────────────────────────────────────────


	%                      End of the simulation timer                        %
	% ─────────────────────────────────────────────────────────────────────────
	simulation_time = toc;
	fprintf('EMT-type simulation completed in:          %.4f seconds\n', simulation_time);
	% ─────────────────────────────────────────────────────────────────────────


	%                                  Plots                                  %
	% ─────────────────────────────────────────────────────────────────────────
	% Sparsity pattern of matrix Gm for small cases
	if n_pi <= 5 && Np <= 5
		figure;
		spy(Gm,'k.',40)
		xlabel('Columns'); ylabel('Rows');
		set(gcf,'pos',[200 200 800 800]); grid on;
	end
	% ─────────────────────────────────────────────────────────────────────────

end