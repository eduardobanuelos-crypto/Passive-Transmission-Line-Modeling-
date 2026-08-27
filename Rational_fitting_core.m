function [Values_ab0,Fs_fit, Fs_fit2] = Rational_fitting_core(Fitting_data, Data, opts)

% This function prepares the impedance functions for their subsequent
% synthesis through a Foster-type RL network using VF-PR, obtaining the
% parameters required for the time-domain simulation of the MTL.

% Inputs:
%     Fitting_data - Structure with the parameters used in the rational fitting:
%         s  = Fitting_data.s;   % Complex frequency
%         Ll = Fitting_data.Ll;  % Transmission-line length
%         Ns = Fitting_data.Ns;  % Number of frequency samples

%     Data  - Structure with the transmission-line parameters in the modal domain:
%         Z_ab0 = Data.Z_ab0; % Modal impedance matrix
%         Y_ab0 = Data.Y_ab0; % Modal admittance matrix
%         Rdc   = Data.Rdc;   % Direct-current resistances

%     opts  - Structure with the options used by NNLS and the plots

% Outputs:
%     Values_ab0 - Structure with the parameters obtained from the rational
%                  fitting for each of the considered modes
%     Fs_fit     - Rational fitting corresponding to the alpha and beta modes
%     Fs_fit2    - Rational fitting corresponding to the zero mode


	%                              Data extraction                            %
	% ─────────────────────────────────────────────────────────────────────────
	Z_ab0  = Data.Z_ab0;
	Y_ab0  = Data.Y_ab0;
	Rdc    = Data.Rdc;
	s      = Fitting_data.s;  % Complex frequency, rad/s
	Ll     = Fitting_data.Ll; % Transmission-line length
	Ns     = Fitting_data.Ns; % Frequency samples
	% ─────────────────────────────────────────────────────────────────────────


	%               Function preparation for Vector Fitting                  %
	% ─────────────────────────────────────────────────────────────────────────
	n_cnds = size(Z_ab0,1);

	if n_cnds == 3

		% Impedances and admittances of the alpha, beta, and zero modes
		Z_ab = reshape(Z_ab0(1,1,:), [1,Ns]); % zero
		Y_ab = reshape(Y_ab0(1,1,:), [1,Ns]);

		Z_0  = reshape(Z_ab0(3,3,:), [1,Ns]); % alpha and beta
		Y_0  = reshape(Y_ab0(3,3,:), [1,Ns]);

		% Preparation of the impedance functions for rational fitting
		Zp_ab = (Z_ab-Rdc(1,1))./s; % Function prepared for fitting (Zp)
		Zp_0  = (Z_0-Rdc(3,1))./s;  % Function prepared for fitting (Zp)

	elseif n_cnds == 6

		% Impedances and admittances of the alpha, beta, zero, and zero2 modes
		Z_ab = reshape(Z_ab0(1,1,:), [1,Ns]); % alpha and beta
		Y_ab = reshape(Y_ab0(1,1,:), [1,Ns]);

		Z_0  = reshape(Z_ab0(3,3,:), [1,Ns]); % zero
		Y_0  = reshape(Y_ab0(3,3,:), [1,Ns]);

		Z_02  = reshape(Z_ab0(4,4,:), [1,Ns]); % zero2
		Y_02  = reshape(Y_ab0(4,4,:), [1,Ns]);

		% Preparation of the impedance functions for rational fitting
		Zp_ab = (Z_ab - Rdc(1,1))./s;  % Function prepared for fitting
		Zp_0  = (Z_0  - Rdc(3,1))./s;  % Function prepared for fitting
		Zp_02 = (Z_02 - Rdc(4,1))./s;  % Function prepared for fitting

	end
	% ─────────────────────────────────────────────────────────────────────────


	%                      Rational fitting with VF and VF-PR                  %
	% ─────────────────────────────────────────────────────────────────────────
	Ka = 2;              % Fitting type
	Plot_index = 1;      % Variable used to store plot data

	if n_cnds == 3

		% Rational fitting of the alpha and beta modes
		[Values_ab, Fs_fit] = VF_VFPR(Fitting_data, opts, Zp_ab, Ka, Rdc(1,1), Ll, Y_ab, Z_ab, Plot_index);
		Plot_index = 2;

		% Rational fitting of the zero mode
		[Values_0, Fs_fit2]  = VF_VFPR(Fitting_data, opts, Zp_0,  Ka, Rdc(3,1), Ll, Y_0,  Z_0, Plot_index);

		% Grouping the results of all modes
		Values_ab0 = [Values_ab; Values_0];

	elseif n_cnds == 6

		% Rational fitting of the alpha and beta modes
		[Values_ab, Fs_fit] = VF_VFPR(Fitting_data, opts, Zp_ab, Ka, Rdc(4,1), Ll, Y_ab, Z_ab, Plot_index);
		Plot_index = 2;

		% Rational fitting of the zero mode
		[Values_0, Fs_fit2]  = VF_VFPR(Fitting_data, opts, Zp_0,  Ka, Rdc(5,1), Ll, Y_0,  Z_0, Plot_index);
		Plot_index = 3;

		% Rational fitting of the second zero mode
		[Values_02, Fs_fit3] = VF_VFPR(Fitting_data, opts, Zp_02, Ka, Rdc(6,1), Ll, Y_02, Z_02, Plot_index);

		% Grouping the results of all modes
		Values_ab0 = [Values_ab; Values_0; Values_02];
	end
	%─────────────────────────────────────────────────────────────────────────



end
