function [Fitting_data, Parameters] = Parameter_calculation_data()

% This function calculates the parameters of a multiconductor
% transmission line (MTL). The program allows the user to select among
% three predefined case studies or enter the geometric and
% electrical data of a line.

% The implementation considers three-phase and double-circuit
% configurations, with or without ground wires.

% Inputs:
%     The function does not require input arguments.

% Outputs:
%     Fitting_data       - Structure with the frequency-sweep parameters
%                          used in the calculation
%     Fitting_parameters - Structure with the electrical parameters of the
%                          line required for rational fitting


	%                              Fitting data                              %
	% ─────────────────────────────────────────────────────────────────────────
	Fitting_data.Ns    = 500;   % Number of samples
	Fitting_data.f_min = 1e-2;  % Minimum frequency, Hz
	Fitting_data.f_max = 1e5;   % Maximum frequency, Hz
	Fitting_data.f = logspace(log10(Fitting_data.f_min), log10(Fitting_data.f_max), Fitting_data.Ns);
	Fitting_data.w = 2.*pi.*Fitting_data.f; 
	Fitting_data.s = 1i.*Fitting_data.w;
	% ─────────────────────────────────────────────────────────────────────────

	%                           Case-study selection                         %
	% ─────────────────────────────────────────────────────────────────────────
	options = { ...
	    'Example 1.- Delta Three-Phase TL', ...
	    'Example 2.- Horizontal Three-Phase TL', ...
	    'Example 3.- Double-Circuit TL', ...
	    'Other Data'};

	[choice, Plots, ok] = Example_menu_editor(options);
	% ─────────────────────────────────────────────────────────────────────────

	if ~ok
	    error('The configuration was canceled.');
	end

	if choice == 1
	%                   Transmission-line data: Case 1                     %
	% ─────────────────────────────────────────────────────────────────────────
	% Format: [type, x_coord, y_coord, n_bundle_cond, radius, bundle_spacing_GMR, resistivity, Length (m)]
	% Case 1: 100-km three-phase transmission line in delta configuration
	    MTL_data = [1, 0,     24.4, 4, 0.01021, 0.4, 2.827433388230814e-08, 100e3;
	                1, 9.27,  28,   4, 0.01021, 0.4, 2.827433388230814e-08, 100e3;
	                1, 18.54, 24.4, 4, 0.01021, 0.4, 2.827433388230814e-08, 100e3];
	% ─────────────────────────────────────────────────────────────────────────

	elseif choice == 2
	%                   Transmission-line data: Case 2                     %
	% ─────────────────────────────────────────────────────────────────────────
	% Case 2: 60.5-km horizontal three-phase transmission line with ground wires
    % Conductor type: 1 = phase, 2 = ground wire
	    MTL_data = [1, 0,     27.9,  4,  0.0155, 0.6, 2.827433388230814e-08, 60.5e3;
	                1, 17.5,  27.9,  4,  0.0155, 0.6, 2.827433388230814e-08, 60.5e3;
	                1, 35,    27.9,  4,  0.0155, 0.6, 2.827433388230814e-08, 60.5e3;
	                2, 4.3,   41.05, 1,   0.008,  0,  2.827433388230814e-08, 60.5e3;
	                2, 30.7,  41.05, 1,   0.008,  0,  2.827433388230814e-08, 60.5e3];
	% ─────────────────────────────────────────────────────────────────────────

	elseif choice == 3
	%                   Transmission-line data: Case 3                     %
	% ─────────────────────────────────────────────────────────────────────────
	% Case 3: 120-km double-circuit three-phase transmission line with ground wires
	    MTL_data = [1, 0,     25.35,  4, 0.01021, 0.4, 2.827433388230814e-08, 120e3;
	                1, 4.7,   35.36,  4, 0.01021, 0.4, 2.827433388230814e-08, 120e3;
	                1, 9.42,  25.35,  4, 0.01021, 0.4, 2.827433388230814e-08, 120e3;
	                1, 22.82, 25.35,  4, 0.01021, 0.4, 2.827433388230814e-08, 120e3;
	                1, 27.52, 35.36,  4, 0.01021, 0.4, 2.827433388230814e-08, 120e3;
	                1, 32.24, 25.35,  4, 0.01021, 0.4, 2.827433388230814e-08, 120e3;
	                2,  6.8,  44.85,  1,  0.008,   0,  2.827433388230814e-08, 120e3;
	                2, 25.44, 44.85,  1,  0.008,   0,  2.827433388230814e-08, 120e3];
	% ─────────────────────────────────────────────────────────────────────────

	elseif choice == 4
	%                   User-defined transmission-line data                 %
	% ─────────────────────────────────────────────────────────────────────────
	% Format:  [type, x_coord, y_coord, n_bundle_cond, radius, bundle_spacing_GMR, resistivity]
	% Conductor type: 1 = phase, 2 = ground wire

	    MTL_data_default = [1, 0,     24.4, 4, 0.01021, 0.4, 2.827433388230814e-08;
	                        1, 9.27,  28,   4, 0.01021, 0.4, 2.827433388230814e-08;
	                        1, 18.54, 24.4, 4, 0.01021, 0.4, 2.827433388230814e-08];

	    line_length_default = 100000; % m

	    [MTL_data, line_length, Fitting_data, ok] = MTL_data_editor( ...
	        MTL_data_default, line_length_default, Fitting_data);

		if ~ok
	        error('The configuration was canceled.');
		end

		MTL_data(:,8) = line_length;
	% ─────────────────────────────────────────────────────────────────────────

	end

	% Store the selection corresponding to plot display
	Aux.Plots = Plots;


	%                          Parameter calculation                        %
	% ─────────────────────────────────────────────────────────────────────────
	[Parameters] = Three_phase_parameter_calculation(Fitting_data, MTL_data, Aux);

	if choice == 4
	    Parameters.length = line_length; % Number of conductors
	elseif choice == 1 || choice == 2 || choice == 3
	    Parameters.length = MTL_data(1,8);
	end
	% ─────────────────────────────────────────────────────────────────────────

	%                         Saving variables to .mat                       %
	% ─────────────────────────────────────────────────────────────────────────
	save('Parameter_calculation.mat', 'Fitting_data', 'Parameters', 'MTL_data');
	% ─────────────────────────────────────────────────────────────────────────

end