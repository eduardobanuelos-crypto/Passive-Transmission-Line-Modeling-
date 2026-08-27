function [Parameters, Z_Trans, Y_Trans, Lo] = Three_phase_parameter_calculation(Fitting_data,Line_data,Aux)

% This function calculates the parameters of the multiconductor transmission line
% and calculates the R, L, and C elements using VFPR and Foster synthesis,
% where passive elements are guaranteed through the use of VFPR and NNLS.
% Ground wires and their subsequent elimination through Kron reduction are included.

	%                              Data extraction                            %
	% ─────────────────────────────────────────────────────────────────────────
	s      = Fitting_data.s;           % Complex frequency, rad/s
	f      = Fitting_data.f;           % Frequency
	w      = Fitting_data.w;           % Frequency, rad/seg
	Ns     = Fitting_data.Ns;          % Frequency samples
	Plot_option     = Aux.Plots;       % Plots 1.- Yes   2.- No
	conductor_type  = Line_data(:,1);  % 1.- Phase   2.- Ground wire
	x_cord = Line_data(:,2);           % x-coordinate of each conductor
	y_cord = Line_data(:,3);           % y-coordinate of each conductor
	n_bundle = Line_data(:,4);         % Number of conductors in the bundle
	r_cnd    = Line_data(:,5);         % Radius of each conductor
	bundle_spacing  = Line_data(:,6);  % Spacing between bundle conductors
	R_al   = Line_data(:,7);           % Aluminum resistivity, Ohm-m
	Lo     = Line_data(1,8);           % Transmission-line length, m
	n_cnds = size(Line_data, 1);       % Number of conductors

	% Phase and ground-wire indices
	phase_idx  = find(conductor_type == 1);
	ground_idx = find(conductor_type == 2);
	n_phase    = length(phase_idx);
	n_ground   = length(ground_idx);
	% ─────────────────────────────────────────────────────────────────────────

	%                         Fundamental constants                          %
	% ─────────────────────────────────────────────────────────────────────────
	Res_earth = 1000;                  % Ground resistivity, Ohm-m
	Muo       = 4*pi*1E-7;             % Magnetic permeability of free space, H/m
	c_light     = 300000000;           % Approximation of c, m/s
	%c_light    = 299792500;           % c, m/s
	Eo        = 1 / (Muo * c_light^2); % Permittivity of free space, F/m
	% ─────────────────────────────────────────────────────────────────────────

	%                   Bundle GMR calculation                              %
	% ─────────────────────────────────────────────────────────────────────────
	gmr_cnds = zeros(n_cnds, 1);
	for nc = 1:n_cnds
		if n_bundle(nc) == 1 % When there is only one conductor
			gmr_cnds(nc) = r_cnd(nc);
        else                 % When there is more than one conductor
			gmr_cnds(nc) = (sqrt(2) * bundle_spacing(nc)^3 * r_cnd(nc))^(1/n_bundle(nc));
		end
	end
	% ─────────────────────────────────────────────────────────────────────────

	%                           Distance matrices                             %
	% ─────────────────────────────────────────────────────────────────────────
	d_c = zeros(n_cnds);    % Conductor-to-conductor distance (m)
	d_ci = zeros(n_cnds);   % Conductor-to-image distance     (m)

	for k = 1:n_cnds
		for m = 1:n_cnds
			if k == m
				% Self-distance is the GMR
				d_ci(k,k) = gmr_cnds(k);
				% For the conductor-to-image distance, the height is used
				d_c(k,k) = 2 * y_cord(k);
			else
				% Distance between conductors
				d_ci(k,m) = sqrt((x_cord(k)-x_cord(m))^2 + (y_cord(k)-y_cord(m))^2);
				% Distance between a conductor and the image of the other conductor
				d_c(k,m) = sqrt((x_cord(k)-x_cord(m))^2 + (y_cord(k)+y_cord(m))^2);
			end
		end
	end
	% ─────────────────────────────────────────────────────────────────────────


	%                 Conductor internal impedance matrix                   %
	% ─────────────────────────────────────────────────────────────────────────
	Al_resistivity = R_al;        % Resistivity vector
	Z_c = zeros(n_cnds, Ns);      % Impedance per conductor, Ohm/m
	Rdc = zeros(n_cnds,1);        % DC resistance, f = 0 Hz, per conductor
	Z_hf = zeros(n_cnds,Ns);      % HF impedance per conductor, Ohm/m

	for n = 1:n_cnds
		% Complex penetration depth in the conductor
		P_c = sqrt(Al_resistivity(n) ./ (s * Muo));
		% High-frequency impedance
		Z_hf(n, :) = (1/n_bundle(n)) * Al_resistivity(n) ./ (2*pi .* r_cnd(n) .* P_c);
		% DC resistance
		Rdc(n, 1) = (1/n_bundle(n)) * Al_resistivity(n) ./ (pi .* (r_cnd(n)^2));
		% Total impedance
		Z_c(n, :) = sqrt(Rdc(n, 1).^2 + Z_hf(n, :).^2);
	end

	% Internal impedance matrix
	Z_cnd = zeros(n_cnds, n_cnds, Ns);
	for k = 1:n_cnds
		for m = 1:n_cnds
			if k == m
				Z_cnd(k, m, :) = Z_c(k, :);
			end
		end
	end
	% ─────────────────────────────────────────────────────────────────────────


	%                 Conductor geometric impedance matrix                  %
	% ─────────────────────────────────────────────────────────────────────────
	L_geo = zeros(n_cnds, n_cnds);
	for k = 1:n_cnds 
		for m = 1:n_cnds
			if k == m
				% Self impedance
				L_geo(k,m) = (Muo/(2*pi)) * log(2*y_cord(k) / gmr_cnds(k));
			else
				% Mutual impedance
				L_geo(k,m) = (Muo/(2*pi)) * log(d_c(k,m) / d_ci(k,m));
			end
		end
	end

	% Geometric impedance matrix
	Z_geo = zeros(n_cnds, n_cnds, Ns);
	for i = 1:Ns
		Z_geo(:, :, i) = s(i) * L_geo;
	end
	% ─────────────────────────────────────────────────────────────────────────


	%                    Ground-return impedance matrix                     %
	% ─────────────────────────────────────────────────────────────────────────
	Z_ear = zeros(n_cnds, n_cnds, Ns);

	% Complex penetration depth in the ground 
	P_t = sqrt(Res_earth ./ (s * Muo));

	for k = 1:n_cnds
		for m = 1:n_cnds
			if k == m
				% Self impedance
				Z_ear(k,k,:) = s .* (Muo/(2*pi)) .* log(1 + P_t./y_cord(k));
			else
				% Distance between conductor and (image + 2P)
				d_ci_2p = sqrt( (x_cord(k) - x_cord(m)).^2 + (y_cord(k) + y_cord(m) + 2 * P_t).^2);
				% Mutual impedance
				Z_ear(k,m,:) = s .* (Muo./(2*pi)) .* log(d_ci_2p ./ d_c(k,m));
			end
		end
	end
	% ─────────────────────────────────────────────────────────────────────────


	%                         Total impedance matrix                        %
	% ─────────────────────────────────────────────────────────────────────────
	Z_tot = Z_cnd + Z_geo + Z_ear;
	% ─────────────────────────────────────────────────────────────────────────


	%                         Total impedance matrix                        %
	% ─────────────────────────────────────────────────────────────────────────
	P_matrix = zeros(n_cnds, n_cnds);
	
	for k = 1:n_cnds
		for m = 1:n_cnds
			if k == m
				% P_ii
				P_matrix(k,m) = (1/(2*pi*Eo)) * log(2*y_cord(k)/gmr_cnds(k));
			else
				% P_ik
				P_matrix(k,m) = (1/(2*pi*Eo)) * log(d_c(k,m)/d_ci(k,m));
			end
		end
	end
	% ─────────────────────────────────────────────────────────────────────────


	%                           Capacitance matrix                           %
	% ─────────────────────────────────────────────────────────────────────────
	%                      C = P^(-1). Y = jw * C = s * C                      
	Cap = inv(P_matrix); % Capacitance matrix C = P^(-1)
	% ─────────────────────────────────────────────────────────────────────────


	%                           Conductance matrix                           %
	% ─────────────────────────────────────────────────────────────────────────
	G = 0; % Not considered
	% ─────────────────────────────────────────────────────────────────────────


	%                         Shunt admittance matrix                       %
	% ─────────────────────────────────────────────────────────────────────────
	Y_tot = zeros(n_cnds, n_cnds, Ns);
	for i = 1:Ns
	  % Y_tot(:, :, i) = s(i) * Cap;
		Y_tot(:, :, i) = s(i) * (P_matrix \ eye(n_cnds));
	end
	% ─────────────────────────────────────────────────────────────────────────


	%                  Kron reduction: Ground-wire elimination              %
	% ─────────────────────────────────────────────────────────────────────────
	if n_ground > 0
		% Split matrices into phase and ground-wire submatrices
		Z_ff = zeros(n_phase,   n_phase,   Ns); % Phase-Phase
		Z_fg = zeros(n_phase,   n_ground, Ns);  % Phase-Ground wire
		Z_gf = zeros(n_ground, n_phase,   Ns);  % Ground wire-Phase
		Z_gg = zeros(n_ground, n_ground, Ns);   % Ground wire-Ground wire
		
		Y_ff = zeros(n_phase,   n_phase,   Ns); % Phase-Phase
		Y_fg = zeros(n_phase,   n_ground, Ns);  % Phase-Ground wire (will be zero)
		Y_gf = zeros(n_ground, n_phase,   Ns);  % Ground wire-Phase (will be zero)
		Y_gg = zeros(n_ground, n_ground, Ns);   % Ground wire-Ground wire
		
		for i = 1:Ns
			% Extract impedance submatrices
			Z_ff(:,:,i) = Z_tot(phase_idx,   phase_idx,   i);
			Z_fg(:,:,i) = Z_tot(phase_idx,   ground_idx, i);
			Z_gf(:,:,i) = Z_tot(ground_idx, phase_idx,   i);
			Z_gg(:,:,i) = Z_tot(ground_idx, ground_idx, i);
			
			% Extract admittance submatrices
			Y_ff(:,:,i) = Y_tot(phase_idx,   phase_idx,   i);
			Y_fg(:,:,i) = Y_tot(phase_idx,   ground_idx, i);
			Y_gf(:,:,i) = Y_tot(ground_idx, phase_idx,   i);
			Y_gg(:,:,i) = Y_tot(ground_idx, ground_idx, i);
		end
		
		% Kron reduction: Z_phase_reduced = Z_ff - Z_fg * inv(Z_gg) * Z_gf
		Z_phase_reduced = zeros(n_phase, n_phase, Ns);
		for i = 1:Ns
			Z_phase_reduced(:,:,i) = Z_ff(:,:,i) - Z_fg(:,:,i) / Z_gg(:,:,i) * Z_gf(:,:,i);
		end
		
		% Kron reduction: Y_phase_reduced = Y_ff - Y_fg * inv(Y_gg) * Y_gf
		Y_phase_reduced = zeros(n_phase, n_phase, Ns);
		for i = 1:Ns
			Y_phase_reduced(:,:,i) = Y_ff(:,:,i); % - Y_fg(:,:,i) / Y_gg(:,:,i) * Y_gf(:,:,i);
		end
		
		% Update total matrices with the reduced matrices
		Z_tot = Z_phase_reduced;
		Y_tot = Y_phase_reduced;
		
		% Update number of conductors
		n_cnds = n_phase;
	end
	% ─────────────────────────────────────────────────────────────────────────


	%             Transposition of impedance and admittance matrices        %
	% ─────────────────────────────────────────────────────────────────────────
	[Z_Trans, Y_Trans] = Transposition(Z_tot, Y_tot, n_cnds, Ns);
	% ─────────────────────────────────────────────────────────────────────────


	%                 Transformation to alpha, beta, and zero               %
	% ─────────────────────────────────────────────────────────────────────────
	[Z_ab0, Y_ab0] = Clarke_transformation(Z_Trans, Y_Trans, Ns, n_cnds);
	% ─────────────────────────────────────────────────────────────────────────


	%                        Modal R, L, C, G parameters                     %
	% ─────────────────────────────────────────────────────────────────────────
	% Reshape changes the dimensions of w from [1,Ns] to [1,1,Ns]
	w_res = reshape(w, [1,1,Ns]);
	% Resistance in ab0
	R_modal = real(Z_ab0);
	% Inductance in ab0
	L_modal = Z_ab0 ./ w_res;
	% Capacitance in ab0
	C_modal = (Y_ab0 ./ w_res);
	% ─────────────────────────────────────────────────────────────────────────


	%                Structure with data for rational fitting               %
	% ─────────────────────────────────────────────────────────────────────────
	Parameters = struct('Z_ab0', Z_ab0, 'Y_ab0', Y_ab0, 'Rdc', Rdc, 'Z_Trans', Z_Trans, 'Y_Trans', Y_Trans, 'Z_tot', Z_tot, 'length', Lo);
	% ─────────────────────────────────────────────────────────────────────────


	%                                  Plots                                 %
	% ─────────────────────────────────────────────────────────────────────────
	if Plot_option == 1 % Plot generation

		% R, L, C abc and ab0 p.u.l.
		figure;
		subplot(1,2,1);
		for k = 1:n_cnds
			for m = 1:n_cnds
				aux1 = abs(squeeze(real(Z_tot(k, m, :))));
				aux2 = abs(squeeze(imag(Z_tot(k, m, :)./w_res )));
				aux3 = abs(squeeze(imag(Y_tot(k, m, :)./w_res )));
				loglog(f,aux1,'r'); hold on;
				loglog(f,aux2,'b');
				loglog(f,aux3,'Color',[0 0.4 0]);
				grid on;
			end
		end
		xlabel('Frequency [Hz]'); ylabel('Magnitude [p.u.l]'); legend('R', 'L', 'C','Location','northwest','Orientation','horizontal');
		xlim([min(f) max(f)]);
		subplot(1,2,2);
		for k = 1:n_cnds
			for m = 1:n_cnds
				aux1 = abs(squeeze(R_modal(k, m, :)));
				aux2 = abs(squeeze(imag(L_modal(k, m, :))));
				aux3 = abs(squeeze(imag(C_modal(k, m, :))));
				if min(aux1) >= 1e-15
					loglog(f,aux1,'r'); hold on;
					loglog(f,aux2,'b');
					loglog(f,aux3,'Color',[0 0.4 0]);
				end
				grid on;
			end
		end
		xlabel('Frequency [Hz]'); ylabel('Magnitude [p.u.l]'); legend('R_{\alpha\beta0}', 'L_{\alpha\beta0}', 'C_{\alpha\beta0}','Location','northwest','Orientation','horizontal');
		xlim([min(f) max(f)]);

		% R, L, C abc and ab0 per Km
		figure;
		subplot(1,2,1);
		for k = 1:n_cnds
			for m = 1:n_cnds
				aux1 = abs(squeeze(real(Z_tot(k, m, :)))) * 1000;
				aux2 = abs(squeeze(imag(Z_tot(k, m, :)./w_res ))) * 1000;
				aux3 = abs(squeeze(imag(Y_tot(k, m, :)./w_res ))) * 1000;
				loglog(f,aux1,'r'); hold on;
				loglog(f,aux2,'b');
				loglog(f,aux3,'Color',[0 0.4 0]); 
				grid on;
			end
		end
		xlabel('Frequency [Hz]'); ylabel('Magnitude'); legend('R [\Omega/Km]', 'L [H/Km]', 'C [F/Km]','Location','best','Orientation','horizontal','fontsize',16);
		xlim([min(f) max(f)]);
		subplot(1,2,2);
		for k = 1:n_cnds
			for m = 1:n_cnds
				aux1 = abs(squeeze(R_modal(k, m, :))) * 1000;
				aux2 = abs(squeeze(imag(L_modal(k, m, :)))) * 1000;
				aux3 = abs(squeeze(imag(C_modal(k, m, :)))) * 1000;
				if min(aux1) >= 1e-15
					loglog(f,aux1,'r'); hold on;
					loglog(f,aux2,'b');
					loglog(f,aux3,'Color',[0 0.4 0]);
				end
				grid on;
			end
		end
		xlabel('Frequency [Hz]'); ylabel('Magnitude'); legend('R_{\alpha\beta0} [\Omega/Km]', 'L_{\alpha\beta0} [H/Km]', 'C_{\alpha\beta0} [F/Km]','Location','best','Orientation','horizontal','fontsize',16);
		xlim([min(f) max(f)]);

	end
	% ─────────────────────────────────────────────────────────────────────────
	
end