function [Data, ok] = Time_simulation_editor(Data_default)

% This function generates a graphical user interface to configure the
% parameters used in the time-domain simulation of the transmission line.

% Inputs:
%     Data_default - Structure with the initial configuration values
%                    for the time-domain simulation

% Outputs:
%     Data - Structure with the simulation parameters selected by
%            the user
%     ok   - Configuration acceptance indicator:
%            true.- Accept    false.- Cancel


	%                         Variable initialization                       %
	% ─────────────────────────────────────────────────────────────────────────
	ok = false;

	% Default values for the solution options
	if ~isfield(Data_default,'Solution_method')
		Data_default.Solution_method = 1;     % 1.- All  2.- SS  3.- EMT  4.- NLT
	end
	if ~isfield(Data_default,'SS_solution_type')
		Data_default.SS_solution_type = 2;    % 1.- Normal  2.- Optimized (sparse)
	end
	if ~isfield(Data_default,'EMT_solution_type')
		Data_default.EMTP_solution_type = 1;  % 1.- Thomas  2.- inv(Gm)*H
	end
	if ~isfield(Data_default,'NLT_window')
		Data_default.NLT_window = 3;          % 1.- Hanning  2.- Wilcox  3.- Blackman  4.- Riesz  5.- Lanczos
	end
	if ~isfield(Data_default,'NLT_damping')
		Data_default.NLT_damping = 1;         % 1.- c=2*dw  2.- c=log(N^2)/T
	end
	if ~isfield(Data_default,'passivity_test')
		Data_default.passivity_test = 2;      % 1.- Yes  2.- No
	end

	Data = Data_default;
	% ─────────────────────────────────────────────────────────────────────────


	%                            Window creation                             %
	% ─────────────────────────────────────────────────────────────────────────
	fig = uifigure('Name', 'Time-Domain Simulation Configuration', ...
				   'Position', [220 80 1050 650]);
	% ─────────────────────────────────────────────────────────────────────────


	%                                  Title                                 %
	% ─────────────────────────────────────────────────────────────────────────
	uilabel(fig, 'Text', 'Time-domain simulation configuration', ...
		'Position', [40 600 780 25], 'FontWeight', 'bold');
	% ─────────────────────────────────────────────────────────────────────────


	%                            Simulation data                            %
	% ─────────────────────────────────────────────────────────────────────────
	uilabel(fig, 'Text', 'Simulation data', 'Position', [40 555 250 25], ...
		'FontWeight', 'bold');

	uilabel(fig, 'Text', 'Time step dt (s):', 'Position', [60 520 210 25]);
	dt_field = uieditfield(fig, 'numeric', ...
		'Value', Data_default.Time_data.dt, 'Position', [280 520 130 25]);

	uilabel(fig, 'Text', 'Final time t_final (s):', 'Position', [60 485 210 25]);
	t_final_field = uieditfield(fig, 'numeric', ...
		'Value', Data_default.Time_data.t_final, 'Position', [280 485 130 25]);

	uilabel(fig, 'Text', 'Final NLT time t_final_NLT (s):', 'Position', [60 450 210 25]);
	t_final_NLT_field = uieditfield(fig, 'numeric', ...
		'Value', Data_default.Time_data.t_final_NLT, 'Position', [280 450 130 25]);
	% ─────────────────────────────────────────────────────────────────────────


	%                       Source and load parameters                      %
	% ─────────────────────────────────────────────────────────────────────────
	uilabel(fig, 'Text', 'Input source and load', 'Position', [40 405 250 25], ...
		'FontWeight', 'bold');

	uilabel(fig, 'Text', 'R_in (Ohm):', 'Position', [60 370 210 25]);
	R_in_field = uieditfield(fig, 'numeric', 'Value', Data_default.Source_values.R_in, ...
		'Position', [280 370 130 25]);

	uilabel(fig, 'Text', 'L_in (H):', 'Position', [60 335 210 25]);
	L_in_field = uieditfield(fig, 'numeric', 'Value', Data_default.Source_values.L_in, ...
		'Position', [280 335 130 25]);

	uilabel(fig, 'Text', 'R_L (Ohm):', 'Position', [60 300 210 25]);
	R_L_field = uieditfield(fig, 'numeric', 'Value', Data_default.Source_values.R_L, ...
		'Position', [280 300 130 25]);

	% ─────────────────────────────────────────────────────────────────────────


	%                               Source type                              %
	% ─────────────────────────────────────────────────────────────────────────
	uilabel(fig, 'Text', 'Source type', 'Position', [460 555 200 25], ...
		'FontWeight', 'bold');

	uilabel(fig, 'Text', 'Source_type:', 'Position', [480 520 100 25]);
	source_type_field = uidropdown(fig, 'Items', {'1.- V_ac', '2.- V_dc'}, ...
		'ItemsData', [1 2], 'Value', Data_default.Source_type, ...
		'Position', [560 520 120 25]);
	% ─────────────────────────────────────────────────────────────────────────


	%                               Source_type AC                                 %
	% ─────────────────────────────────────────────────────────────────────────
	uilabel(fig, 'Text', 'Source_type AC', 'Position', [460 475 200 25], 'FontWeight', 'bold');
	uilabel(fig, 'Text', 'Vac_a:', 'Position', [480 440 80 25]);
	Vac_a_field = uieditfield(fig, 'numeric', 'Value', Data_default.Voltage_values.Vac(1,1), 'Position', [560 440 100 25]);
	uilabel(fig, 'Text', 'Vac_b:', 'Position', [480 405 80 25]);
	Vac_b_field = uieditfield(fig, 'numeric', 'Value', Data_default.Voltage_values.Vac(2,1), 'Position', [560 405 100 25]);
	uilabel(fig, 'Text', 'Vac_c:', 'Position', [480 370 80 25]);
	Vac_c_field = uieditfield(fig, 'numeric', 'Value', Data_default.Voltage_values.Vac(3,1), 'Position', [560 370 100 25]);
	uilabel(fig, 'Text', 'Ang_a (°):', 'Position', [480 330 80 25]);
	angle_a_field = uieditfield(fig,   'numeric', 'Value', Data_default.Voltage_values.Angles(1,1), 'Position', [560 330 100 25]);
	uilabel(fig, 'Text', 'Ang_b (°):', 'Position', [480 295 80 25]);
	angle_b_field = uieditfield(fig,   'numeric', 'Value', Data_default.Voltage_values.Angles(2,1), 'Position', [560 295 100 25]);
	uilabel(fig, 'Text', 'Ang_c (°):', 'Position', [480 260 80 25]);
	angle_c_field = uieditfield(fig,   'numeric', 'Value', Data_default.Voltage_values.Angles(3,1), 'Position', [560 260 100 25]);
	uilabel(fig, 'Text', 'General_angle (°):', 'Position', [480 220 80 25]);
	general_angle_field = uieditfield(fig, 'numeric', 'Value', Data_default.Voltage_values.General_angle, 'Position', [560 220 100 25]);
	% ─────────────────────────────────────────────────────────────────────────


	%                               Source_type DC                                 %
	% ─────────────────────────────────────────────────────────────────────────
	uilabel(fig, 'Text', 'Source_type DC', 'Position', [40 215 250 25], 'FontWeight', 'bold');
	uilabel(fig, 'Text', 'Vdc_a:', 'Position', [60 180 80 25]);
	Vdc_a_field = uieditfield(fig, 'numeric', 'Value', Data_default.Voltage_values.Vdc(1,1), 'Position', [140 180 100 25]);
	uilabel(fig, 'Text', 'Vdc_b:', 'Position', [60 145 80 25]);
	Vdc_b_field = uieditfield(fig, 'numeric', 'Value', Data_default.Voltage_values.Vdc(2,1), 'Position', [140 145 100 25]);
	uilabel(fig, 'Text', 'Vdc_c:', 'Position', [60 110 80 25]);
	Vdc_c_field = uieditfield(fig, 'numeric', 'Value', Data_default.Voltage_values.Vdc(3,1), 'Position', [140 110 100 25]);
	% ─────────────────────────────────────────────────────────────────────────


	%                            Solution methods                           %
	% ─────────────────────────────────────────────────────────────────────────
	uilabel(fig, 'Text', 'Solution methods', 'Position', [730 555 250 25], ...
		'FontWeight', 'bold');

	uilabel(fig, 'Text', 'Run:', 'Position', [750 520 100 25]);
	solution_method_field = uidropdown(fig, ...
		'Items', {'All', 'State Space', 'EMT', 'NLT'}, ...
		'ItemsData', [1 2 3 4], 'Value', Data_default.Solution_method, ...
		'Position', [840 520 165 25]);

	% State-space representation options
	uilabel(fig, 'Text', 'State Space', 'Position', [730 475 250 25], 'FontWeight', 'bold');
	uilabel(fig, 'Text', 'Method:', 'Position', [750 440 80 25]);
	SS_method_field = uidropdown(fig, ...
		'Items', {'Normal', 'Optimized (sparse)'}, 'ItemsData', [1 2], ...
		'Value', Data_default.SS_solution_type, 'Position', [840 440 165 25]);

	% EMT-type solution options
	uilabel(fig, 'Text', 'EMT-type solution', 'Position', [730 395 250 25], 'FontWeight', 'bold');
	uilabel(fig, 'Text', 'Method:', 'Position', [750 360 80 25]);
	EMT_method_field = uidropdown(fig, ...
		'Items', {'Thomas Algorithm', 'Inverse of G_m'}, 'ItemsData', [1 2], ...
		'Value', Data_default.EMTP_solution_type, 'Position', [840 360 165 25]);

	% Numerical Laplace Transform options
	uilabel(fig, 'Text', 'Numerical Laplace Transform', 'Position', [730 315 280 25], 'FontWeight', 'bold');
	uilabel(fig, 'Text', 'Window:', 'Position', [750 280 80 25]);
	window_field = uidropdown(fig, ...
		'Items', {'Hanning', 'Wilcox', 'Blackman', 'Riesz', 'Lanczos'}, ...
		'ItemsData', [1 2 3 4 5], 'Value', Data_default.NLT_window, ...
		'Position', [840 280 165 25]);

	uilabel(fig, 'Text', 'Damping:', 'Position', [730 240 110 25]);
	damping_field = uidropdown(fig, 'Items', {'c = 2*dw', 'c = log(N^2)/T'}, ...
		'ItemsData', [1 2], 'Value', Data_default.NLT_damping, ...
		'Position', [840 240 165 25]);
	% ─────────────────────────────────────────────────────────────────────────


	%                             Passivity test                            %
	% ─────────────────────────────────────────────────────────────────────────
	uilabel(fig, 'Text', 'Passivity test', 'Position', [730 195 250 25], ...
		'FontWeight', 'bold');

	uilabel(fig, 'Text', 'Run test?', 'Position', [730 160 100 25]);
	passivity_field = uidropdown(fig, 'Items', {'Yes', 'No'}, ...
		'ItemsData', [1 2], 'Value', Data_default.passivity_test, ...
		'Position', [840 160 165 25]);

	% Enable only the options associated with the selected method
	solution_method_field.ValueChangedFcn = @(src,event) update_methods( ...
		solution_method_field, SS_method_field, EMT_method_field, window_field, damping_field);
	update_methods(solution_method_field, SS_method_field, EMT_method_field, window_field, damping_field);
	% ─────────────────────────────────────────────────────────────────────────


	%                                 Buttons                                %
	% ─────────────────────────────────────────────────────────────────────────
	uibutton(fig, 'push', 'Text', 'Accept', 'Position', [790 75 100 35], ...
		'ButtonPushedFcn', @(btn, event) accept_simulation_data( ...
			fig, Data_default, dt_field, t_final_field, t_final_NLT_field, ...
			R_in_field, L_in_field, R_L_field, source_type_field, ...
			Vac_a_field, Vac_b_field, Vac_c_field, ...
			angle_a_field, angle_b_field, angle_c_field, general_angle_field, ...
			Vdc_a_field, Vdc_b_field, Vdc_c_field, ...
			solution_method_field, SS_method_field, EMT_method_field, window_field, damping_field, ...
			passivity_field));

	uibutton(fig, 'push', 'Text', 'Cancel', 'Position', [910 75 100 35], ...
		'ButtonPushedFcn', @(btn, event) cancel_simulation(fig));

	% Prevent uncontrolled window closing
	fig.CloseRequestFcn = @(src, event) cancel_simulation(fig);
	% ─────────────────────────────────────────────────────────────────────────


	%                       Wait and data retrieval                         %
	% ─────────────────────────────────────────────────────────────────────────
	% Wait for the user to press Accept or Cancel
	uiwait(fig);

	% Retrieve data
	if isvalid(fig)
		ok = getappdata(fig, 'ok');
		if ok
			Data = getappdata(fig, 'Data');
		end
		delete(fig);
	end
	% ─────────────────────────────────────────────────────────────────────────
end


function accept_simulation_data(fig, Data_default, ...
	dt_field, t_final_field, t_final_NLT_field, ...
	R_in_field, L_in_field, R_L_field, source_type_field, ...
	Vac_a_field, Vac_b_field, Vac_c_field, ...
	angle_a_field, angle_b_field, angle_c_field, general_angle_field, ...
	Vdc_a_field, Vdc_b_field, Vdc_c_field, ...
	solution_method_field, SS_method_field, EMT_method_field, window_field, damping_field, ...
	passivity_field)

	%                          Parameter reading                            %
	% ─────────────────────────────────────────────────────────────────────────
	% Read simulation data
	dt          = dt_field.Value;
	t_final     = t_final_field.Value;
	t_final_NLT = t_final_NLT_field.Value;

	% Read source and load parameters
	R_in = R_in_field.Value;
	L_in = L_in_field.Value;
	R_L  = R_L_field.Value;

	% Read source type
	Source_type = source_type_field.Value;

	% Read AC source
	Vac = [Vac_a_field.Value;
		   Vac_b_field.Value;
		   Vac_c_field.Value];

	Angles = [angle_a_field.Value;
		      angle_b_field.Value;
		      angle_c_field.Value];

	General_angle = general_angle_field.Value;

	% Read DC source
	Vdc = [Vdc_a_field.Value;
		   Vdc_b_field.Value;
		   Vdc_c_field.Value];

	% Read solution methods
	Solution_method    = solution_method_field.Value;
	SS_solution_type   = SS_method_field.Value;
	EMT_solution_type  = EMT_method_field.Value;
	NLT_window         = window_field.Value;
	NLT_damping        = damping_field.Value;
	passivity_test     = passivity_field.Value;
	% ─────────────────────────────────────────────────────────────────────────


	%                              Validations                               %
	% ─────────────────────────────────────────────────────────────────────────
	if dt <= 0
		uialert(fig, 'The time step dt must be greater than zero.', 'Error');
		return;
	end
	if t_final <= 0
		uialert(fig, 'The final time t_final must be greater than zero.', 'Error');
		return;
	end
	if t_final_NLT <= 0
		uialert(fig, 'The final time t_final_NLT must be greater than zero.', 'Error');
		return;
	end
	if t_final <= dt
		uialert(fig, 'The final time t_final must be greater than dt.', 'Error');
		return;
	end
	if t_final_NLT <= dt
		uialert(fig, 'The final time t_final_NLT must be greater than dt.', 'Error');
		return;
	end
	if R_in < 0
		uialert(fig, 'The resistance R_in cannot be negative.', 'Error');
		return;
	end
	if L_in < 0
		uialert(fig, 'The inductance L_in cannot be negative.', 'Error');
		return;
	end
	if R_L <= 0
		uialert(fig, 'The load resistance R_L must be greater than zero.', 'Error');
		return;
	end
	% ─────────────────────────────────────────────────────────────────────────


	%                               Data storage                              %
	% ─────────────────────────────────────────────────────────────────────────
	Data = Data_default;

	% Simulation
	Data.Time_data.dt          = dt;
	Data.Time_data.t_final     = t_final;
	Data.Time_data.t_final_NLT = t_final_NLT;

	Data.Time_data.t  = 0:dt:t_final-dt;
	Data.Time_data.N  = length(Data.Time_data.t);
	Data.Time_data.t2 = 0:dt:t_final_NLT-dt;
	Data.Time_data.N2 = length(Data.Time_data.t2);

	% Source and load
	Data.Source_values.R_in = R_in;
	Data.Source_values.L_in = L_in;
	Data.Source_values.R_L  = R_L;

	% Source type
	Data.Source_type = Source_type;

	% Source_type AC
	Data.Voltage_values.Vac    = Vac;
	Data.Voltage_values.Angles = Angles;
	Data.Voltage_values.General_angle = General_angle;

	% Source_type DC
	Data.Voltage_values.Vdc = Vdc;

	% Solution methods
	Data.Solution_method    = Solution_method;
	Data.SS_solution_type   = SS_solution_type;
	Data.EMT_solution_type  = EMT_solution_type;
	Data.NLT_window         = NLT_window;
	Data.NLT_damping        = NLT_damping;
	Data.passivity_test     = passivity_test;

	setappdata(fig, 'Data', Data);
	setappdata(fig, 'ok', true);
	uiresume(fig);
	% ─────────────────────────────────────────────────────────────────────────
end


function update_methods(solution_method_field, SS_method_field, EMT_method_field, window_field, damping_field)

	%                      Enabling solution methods                       %
	% ─────────────────────────────────────────────────────────────────────────
	% Initially disable all specific options
	SS_method_field.Enable  = 'off';
	EMT_method_field.Enable = 'off';
	window_field.Enable  = 'off';
	damping_field.Enable = 'off';

	% Enable only the required options
	if solution_method_field.Value == 1       % All
		SS_method_field.Enable     = 'on';
		EMT_method_field.Enable   = 'on';
		window_field.Enable  = 'on';
		damping_field.Enable = 'on';
	elseif solution_method_field.Value == 2   % State Space
		SS_method_field.Enable      = 'on';
	elseif solution_method_field.Value == 3   % EMTP
		EMT_method_field.Enable    = 'on';
	elseif solution_method_field.Value == 4   % NLT
		window_field.Enable  = 'on';
		damping_field.Enable = 'on';
	end
	% ─────────────────────────────────────────────────────────────────────────
end


function cancel_simulation(fig)

	%                       Configuration cancellation                      %
	% ─────────────────────────────────────────────────────────────────────────
	setappdata(fig, 'ok', false);
	uiresume(fig);
	% ─────────────────────────────────────────────────────────────────────────
end
