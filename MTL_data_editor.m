function [MTL_data, line_length, Fitting_data, ok] = MTL_data_editor(Data_default, line_length_default, Fitting_data_default)

% This function generates a graphical user interface to configure the
% geometric and electrical data of a multiconductor transmission
% line, as well as its length and the parameters of the frequency
% sweep subsequently used in the parameter calculation.

% Inputs:
%     Data_default            - Matrix with the initial conductor data
%     line_length_default     - Initial transmission-line length, m
%     Fitting_data_default    - Structure with the initial parameters of the
%                               frequency sweep:
%     Ns    = Fitting_data_default.Ns;    % Number of samples
%     f_min = Fitting_data_default.f_min; % Minimum frequency, Hz
%     f_max = Fitting_data_default.f_max; % Maximum frequency, Hz

% Outputs:
%     MTL_data     - Matrix with the configured conductor data
%     line_length  - Transmission-line length, m
%     Fitting_data - Structure with the frequency vector and parameters
%                    associated with the frequency sweep
%     ok           - Data acceptance indicator:
%                    true.- Accept    false.- Cancel

	%                       Variable initialization                      %
	% ─────────────────────────────────────────────────────────────────────
	ok = false;
	MTL_data = Data_default;
	line_length = line_length_default;
	Fitting_data = Fitting_data_default;
	% ─────────────────────────────────────────────────────────────────────

	
	%                           Window creation                           %
	% ─────────────────────────────────────────────────────────────────────
	% Create window
	fig = uifigure('Name', 'Multiconductor Transmission Line Configuration', ...
				   'Position', [300 150 1000 520]);

	% Variable used to store selected rows
	setappdata(fig, 'SelectedRows', []);
	% ─────────────────────────────────────────────────────────────────────


	%                         General line data                          %
	% ─────────────────────────────────────────────────────────────────────

	uilabel(fig, 'Text', 'Line length (m):', 'Position', [30 465 160 25]);
	line_length_field = uieditfield(fig, 'numeric', 'Value', line_length_default, 'Position', [160 465 100 25]);
	% ─────────────────────────────────────────────────────────────────────


	%                       Frequency-fitting data                       %
	% ─────────────────────────────────────────────────────────────────────
	uilabel(fig, 'Text', 'Number of samples Ns:', 'Position', [330 465 160 25]);

	Ns_field = uieditfield(fig, 'numeric', 'Value', Fitting_data_default.Ns, ...
		'RoundFractionalValues', 'on', ...
		'Position', [480 465 100 25]);

	uilabel(fig, 'Text', 'f min (Hz):', 'Position', [630 465 80 25]);

	f_min_field = uieditfield(fig, 'numeric', 'Value', Fitting_data_default.f_min, ...
		'Position', [700 465 90 25]);

	uilabel(fig, 'Text', 'f max (Hz):', 'Position', [820 465 80 25]);

	f_max_field = uieditfield(fig, 'numeric', 'Value', Fitting_data_default.f_max, ...
		'Position', [890 465 80 25]);
	% ─────────────────────────────────────────────────────────────────────


	%                          Conductor table                           %
	% ─────────────────────────────────────────────────────────────────────
	column_names = {'Type', 'x (m)', 'y (m)', ...
						'No. bundle cond.', 'Radius (m)', 'Bundle dist./GMR (m)', 'Resistivity (Ohm-m)'};

	conductor_table = uitable(fig, 'Data', Data_default, ...
		'ColumnName', column_names, ...
		'ColumnEditable', true(1, 7), ...
		'ColumnFormat', {'numeric', 'numeric', 'numeric', ...
						 'numeric', 'numeric', 'numeric', 'numeric'}, ...
		'Position', [30 130 940 300]);

	% Store selected row
	conductor_table.CellSelectionCallback = @(src, event) select_row(fig, event);

	% Instructions
	uilabel(fig, 'Text', 'Type: 1 = phase, 2 = ground wire. ', ...
		'Position', [30 435 600 25]);
	% ─────────────────────────────────────────────────────────────────────


	%                               Buttons                              %
	% ─────────────────────────────────────────────────────────────────────
	uibutton(fig, 'push', 'Text', 'Add phase', 'Position', [30 70 120 30], ...
		'ButtonPushedFcn', @(btn, event) add_phase(conductor_table));

	uibutton(fig, 'push', 'Text', 'Add ground wire', ...
		'Position', [165 70 120 30], ...
		'ButtonPushedFcn', @(btn, event) add_ground_wire(conductor_table));

	uibutton(fig, 'push', ...
		'Text', 'Delete selected', ...
		'Position', [300 70 160 30], ...
		'ButtonPushedFcn', @(btn, event) delete_row(fig, conductor_table));

	uibutton(fig, 'push', ...
		'Text', 'Accept', ...
		'Position', [725 70 110 30], ...
		'ButtonPushedFcn', @(btn, event) accept_data( ...
			fig, conductor_table, line_length_field, Ns_field, f_min_field, f_max_field));

	uibutton(fig, 'push', ...
		'Text', 'Cancel', ...
		'Position', [860 70 110 30], ...
		'ButtonPushedFcn', @(btn, event) cancel_editor(fig));

	% Prevent uncontrolled window closing
	fig.CloseRequestFcn = @(src, event) cancel_editor(fig);
	% ─────────────────────────────────────────────────────────────────────


	%                         Data retrieval                             %
	% ─────────────────────────────────────────────────────────────────────
	% Wait for the user to press Accept or Cancel
	uiwait(fig);

	% Retrieve data
	if isvalid(fig)
		ok = getappdata(fig, 'ok');

		if ok
			MTL_data    = getappdata(fig, 'MTL_data');
			line_length     = getappdata(fig, 'line_length');
			Fitting_data = getappdata(fig, 'Fitting_data');
		end

		delete(fig);
	end
	% ─────────────────────────────────────────────────────────────────────

end


function select_row(fig, event)

	%                         Table-row selection                        %
	% ─────────────────────────────────────────────────────────────────────
	if isempty(event.Indices)
		setappdata(fig, 'SelectedRows', []);
	else
		rows = unique(event.Indices(:, 1));
		setappdata(fig, 'SelectedRows', rows);
	end
	% ─────────────────────────────────────────────────────────────────────

end


function add_phase(conductor_table)

	%                         Add phase conductor                        %
	% ─────────────────────────────────────────────────────────────────────
	Data = conductor_table.Data;

	% Default row for a phase conductor
	new_phase = [1, 0, 24, 4, 0.01021, 0.4, 2.827433388230814e-08];
	conductor_table.Data = [Data; new_phase];
	% ─────────────────────────────────────────────────────────────────────

end


function add_ground_wire(conductor_table)

	%                         Add ground wire                            %
	% ─────────────────────────────────────────────────────────────────────
	Data = conductor_table.Data;

	% Default row for a ground wire
	new_ground_wire = [2, 4, 40, 1, 0.008, 0, 2.827433388230814e-08];

	conductor_table.Data = [Data; new_ground_wire];
	% ─────────────────────────────────────────────────────────────────────

end


function delete_row(fig, conductor_table)

	%                           Row deletion                              %
	% ─────────────────────────────────────────────────────────────────────
	Data = conductor_table.Data;
	rows = getappdata(fig, 'SelectedRows');

	if isempty(rows)
		uialert(fig, ...
			'First select one or more rows in the table.', ...
			'No selection');
		return;
	end

	Data(rows, :) = [];

	if isempty(Data)
		uialert(fig, 'The table cannot be left empty.', 'Insufficient data');
		return;
	end

	conductor_table.Data = Data;
	setappdata(fig, 'SelectedRows', []);
	% ─────────────────────────────────────────────────────────────────────

end


function accept_data(fig, conductor_table, line_length_field, Ns_field, f_min_field, f_max_field)

	%                        Reading entered data                        %
	% ─────────────────────────────────────────────────────────────────────
	Data = conductor_table.Data;
	line_length = line_length_field.Value;
	Ns       = Ns_field.Value;
	f_min    = f_min_field.Value;
	f_max    = f_max_field.Value;
	% ─────────────────────────────────────────────────────────────────────


	%                           Table validations                         %
	% ─────────────────────────────────────────────────────────────────────

	if isempty(Data)
		uialert(fig, 'The table cannot be empty.', 'Error');
		return;
	end

	if size(Data, 2) ~= 7
		uialert(fig, ...
			'The matrix must have exactly 7 columns.', ...
			'Error');
		return;
	end

	if any(isnan(Data), 'all')
		uialert(fig, ...
			'There are empty or nonnumeric values in the table.', ...
			'Error');
		return;
	end

	if any(Data(:, 1) ~= 1 & Data(:, 1) ~= 2)
		uialert(fig, ...
			'The first column can only contain type 1 = phase or type 2 = ground wire.', ...
			'Error');
		return;
	end

	if any(Data(:, 4) < 1)
		uialert(fig, ...
			'The number of conductors in the bundle must be greater than or equal to 1.', ...
			'Error');
		return;
	end

	if any(Data(:, 5) <= 0)
		uialert(fig, ...
			'The conductor radius must be greater than zero.', ...
			'Error');
		return;
	end

	% Check conductor resistivity
	if any(Data(:, 7) <= 0)
		uialert(fig, ...
			'The conductor resistivity must be greater than zero.', ...
			'Error');
		return;
	end
	% ─────────────────────────────────────────────────────────────────────


	%                  Length and frequency-range validations             %
	% ─────────────────────────────────────────────────────────────────────
	if line_length <= 0
		uialert(fig, ...
			'The line length must be greater than zero.', ...
			'Error');
		return;
	end

	if Ns < 2
		uialert(fig, ...
			'The number of samples Ns must be greater than or equal to 2.', ...
			'Error');
		return;
	end

	if mod(Ns, 1) ~= 0
		uialert(fig, ...
			'The number of samples Ns must be an integer.', ...
			'Error');
		return;
	end

	if f_min <= 0
		uialert(fig, ...
			'The minimum frequency must be greater than zero.', ...
			'Error');
		return;
	end

	if f_max <= 0
		uialert(fig, ...
			'The maximum frequency must be greater than zero.', ...
			'Error');
		return;
	end

	if f_max <= f_min
		uialert(fig, ...
			'The maximum frequency must be greater than the minimum frequency.', ...
			'Error');
		return;
	end
	% ─────────────────────────────────────────────────────────────────────


	%                     Frequency-vector construction                   %
	% ─────────────────────────────────────────────────────────────────────
	Fitting_data.Ns    = Ns;
	Fitting_data.f_min = f_min;
	Fitting_data.f_max = f_max;

	Fitting_data.f = logspace(log10(f_min), log10(f_max), Ns);
	Fitting_data.w = 2.*pi.*Fitting_data.f;
	Fitting_data.s = 1i.*Fitting_data.w;
	% ─────────────────────────────────────────────────────────────────────


	%                       Data storage and closing                      %
	% ─────────────────────────────────────────────────────────────────────
	% Store data and close window
	setappdata(fig, 'MTL_data', Data);
	setappdata(fig, 'line_length', line_length);
	setappdata(fig, 'Fitting_data', Fitting_data);
	setappdata(fig, 'ok', true);
	uiresume(fig);
	% ─────────────────────────────────────────────────────────────────────

end


function cancel_editor(fig)

	%                       Configuration cancellation                    %
	% ─────────────────────────────────────────────────────────────────────
	setappdata(fig, 'ok', false);
	uiresume(fig);
	% ─────────────────────────────────────────────────────────────────────
end