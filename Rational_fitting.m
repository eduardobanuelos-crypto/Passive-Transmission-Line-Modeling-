function [Values_ab0, Fitting_data]=Rational_fitting()

% This function performs the rational fitting of the previously
% calculated transmission-line parameters using Vector Fitting with
% Real Poles (VFPR) and the Non-Negative Least Squares (NNLS) algorithm
% selected by the user.

% Inputs:
%     The function does not require input arguments.

% Outputs:
%     Values_ab0    - Structure with the parameters obtained through
%                     rational fitting for the considered modes
%     Fitting_data  - Structure with the parameters used in the fitting


	%                       Reading and variable extraction                  %
	% ─────────────────────────────────────────────────────────────────────────
	load('Parameter_calculation.mat', 'Fitting_data', 'Parameters');

	Fitting_data.Ll  = Parameters.length; % Transmission-line length (m)
	% ─────────────────────────────────────────────────────────────────────────


	%                         Default fitting values                         %
	% ─────────────────────────────────────────────────────────────────────────
	Fitting_data.Np         = 10;  % Number of poles
	Fitting_data.iterations = 15;  % Number of VF and VFPR iterations
	Fitting_data.M          = 100; % Number of cascaded Pi circuits
	opts.Plots              = 1;   % Show plots? 1.- Yes, 2.- No
	opts.NNLS               = 5;   % NNLS type
	% ─────────────────────────────────────────────────────────────────────────


	%                    Rational fitting configuration                     %
	% ─────────────────────────────────────────────────────────────────────────
	[Fitting_data, opts] = VF_fitting_editor(Fitting_data, opts);
	% ─────────────────────────────────────────────────────────────────────────


	%                       Rational fitting execution                       %
	% ─────────────────────────────────────────────────────────────────────────
	[Values_ab0,~,~] = Rational_fitting_core(Fitting_data, Parameters, opts);
	 Values_ab0(1).M = Fitting_data.M; Values_ab0(2).M = Values_ab0(1).M;
	% ─────────────────────────────────────────────────────────────────────────


	%                         Saving variables to .mat                       %
	% ─────────────────────────────────────────────────────────────────────────
	save('Simulation_parameters.mat', 'Values_ab0', 'Fitting_data');
	% ─────────────────────────────────────────────────────────────────────────

end