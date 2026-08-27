function [Fitting_data, opts] = VF_fitting_editor(Fitting_data_default, opts_default)

% This function generates a graphical user interface to configure the
% main parameters used in rational fitting through Vector
% Fitting (VF) and Vector Fitting with Real Poles (VFPR).

% Inputs:
%     Fitting_data_default - Structure with the initial values of the
%                            main fitting parameters:
%         Np  = Fitting_data_default.Np;  % Number of poles
%         iterations = Fitting_data_default.iterations; % Number of VF/VFPR iterations
%         n   = Fitting_data_default.n;   % Number of cascaded Pi circuits

%     opts_default - Structure with the initial fitting options:
%         Plots = opts_default.Plots; % 1.- Yes   2.- No
%         NNLS     = opts_default.NNLS;     % NNLS algorithm type

% Outputs:
%     Fitting_data - Structure with the fitting parameters selected
%                    by the user
%     opts         - Structure with the plot options and NNLS type


	%                           Data initialization                          %
	% ─────────────────────────────────────────────────────────────────────
	Fitting_data = Fitting_data_default;
	opts = opts_default;
	% ─────────────────────────────────────────────────────────────────────


	%                            Window creation                              %
	% ─────────────────────────────────────────────────────────────────────
	fig = uifigure('Name', 'VF and VFPR Fitting Configuration', ...
				   'Position', [400 250 620 390]);
	% ─────────────────────────────────────────────────────────────────────


	%                         Main fitting parameters                         %
	% ─────────────────────────────────────────────────────────────────────
	uilabel(fig, 'Text', 'Number of poles:', 'Position', [40 315 220 25]);

	Np_field = uieditfield(fig, 'numeric', 'Value', Fitting_data_default.Np, ...
		'RoundFractionalValues', 'on', 'Position', [290 315 120 25]);

	uilabel(fig, 'Text', 'Number of VF/VFPR iterations:', 'Position', [40 270 220 25]);

	iterations_field = uieditfield(fig, 'numeric', 'Value', Fitting_data_default.iterations, ...
		'RoundFractionalValues', 'on', 'Position', [290 270 120 25]);

	uilabel(fig, 'Text', 'Number of cascaded Pi circuits:', 'Position', [40 225 230 25]);

	M_field = uieditfield(fig, 'numeric', 'Value', Fitting_data_default.M, ...
		'RoundFractionalValues', 'on', 'Position', [290 225 120 25]);
	% ─────────────────────────────────────────────────────────────────────


	%                              Plot options                              %
	% ─────────────────────────────────────────────────────────────────────
	uilabel(fig, 'Text', 'Show plots:', 'Position', [40 170 220 25]);

	plots_field = uidropdown(fig, 'Items', {'1.- Yes', '2.- No'}, ...
		'ItemsData', [1 2], 'Value', opts_default.Plots, 'Position', [290 170 180 25]);
	% ─────────────────────────────────────────────────────────────────────


	%                              NNLS type                                 %
	% ─────────────────────────────────────────────────────────────────────
	uilabel(fig, 'Text', 'NNLS type:', 'Position', [40 125 220 25]);

	NNLS_field = uidropdown(fig, ...
		'Items', {'1.- NNLS Matlab', ...
				  '2.- FNNLS', ...
				  '3.- TNT-NN', ...
				  '4.- L&W-NNLS', ...
				  '5.- BW-NNLS'}, ...
		'ItemsData', [1 2 3 4 5], ...
		'Value', opts_default.NNLS, ...
		'Position', [290 125 230 25]);
	% ─────────────────────────────────────────────────────────────────────


	%                          Explanatory text                             %
	% ─────────────────────────────────────────────────────────────────────
	uilabel(fig, 'Text', 'Configure the main rational fitting parameters.', ...
		'Position', [40 350 520 25], 'FontWeight', 'bold');
	% ─────────────────────────────────────────────────────────────────────


	%                               Buttons                                  %
	% ─────────────────────────────────────────────────────────────────────
	uibutton(fig, 'push', 'Text', 'Accept', 'Position', [340 45 110 35], ...
		'ButtonPushedFcn', @(btn, event) accept_data( ...
			fig, Fitting_data_default, Np_field, iterations_field, M_field, ...
			plots_field, NNLS_field));

	uibutton(fig, 'push', 'Text', 'Cancel', 'Position', [470 45 110 35], ...
		'ButtonPushedFcn', @(btn, event) cancel_editor(fig));

	% Prevent uncontrolled window closing
	fig.CloseRequestFcn = @(src, event) cancel_editor(fig);
	% ─────────────────────────────────────────────────────────────────────


	%                       Wait and data retrieval                         %
	% ─────────────────────────────────────────────────────────────────────
	% Wait for the user to press Accept or Cancel
	uiwait(fig);

	% Retrieve data
	if isvalid(fig)

		ok = getappdata(fig, 'ok');

		if ok
			Fitting_data = getappdata(fig, 'Fitting_data');
			opts = getappdata(fig, 'opts');
		end

		delete(fig);

	end
	% ─────────────────────────────────────────────────────────────────────

end


%                         Function to accept data                         %
% ─────────────────────────────────────────────────────────────────────────
function accept_data(fig, Fitting_data_default, Np_field, iterations_field, M_field, plots_field, NNLS_field)

	Np  = Np_field.Value;
	iterations = iterations_field.Value;
	M   = M_field.Value;

	Plots = plots_field.Value;
	NNLS     = NNLS_field.Value;

	%                              Validations                               %
	% ─────────────────────────────────────────────────────────────────────
	if Np < 1
		uialert(fig, 'The number of poles must be greater than or equal to 1.', 'Error');
		return;
	end

	if mod(Np, 1) ~= 0
		uialert(fig, 'The number of poles must be an integer.', 'Error');
		return;
	end

	if iterations < 1
		uialert(fig, 'The number of iterations must be greater than or equal to 1.', 'Error');
		return;
	end

	if mod(iterations, 1) ~= 0
		uialert(fig, 'The number of iterations must be an integer.', 'Error');
		return;
	end

	if M < 1
		uialert(fig, 'The number of cascaded Pi circuits must be greater than or equal to 1.', 'Error');
		return;
	end

	if mod(M, 1) ~= 0
		uialert(fig, 'The number of cascaded Pi circuits must be an integer.', 'Error');
		return;
	end
	% ─────────────────────────────────────────────────────────────────────


	%                              Save data                                 %
	% ─────────────────────────────────────────────────────────────────────
	Fitting_data            = Fitting_data_default;
	Fitting_data.Np         = Np;
	Fitting_data.iterations = iterations;
	Fitting_data.M          = M;

	opts.Plots = Plots;
	opts.NNLS  = NNLS;

	setappdata(fig, 'Fitting_data', Fitting_data);
	setappdata(fig, 'opts', opts);
	setappdata(fig, 'ok', true);

	uiresume(fig);
	% ─────────────────────────────────────────────────────────────────────

end
% ─────────────────────────────────────────────────────────────────────────


%                           Function to cancel                           %
% ─────────────────────────────────────────────────────────────────────────
function cancel_editor(fig)
	setappdata(fig, 'ok', false);
	uiresume(fig);
end
% ─────────────────────────────────────────────────────────────────────────
