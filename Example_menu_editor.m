function [choice, Plots, ok] = Example_menu_editor(options)

% This function generates a graphical user interface to select one of
% the examples available for simulation and define whether the associated
% plots will be displayed.

% Inputs:
%     options - Array containing the names of the available examples

% Outputs:
%     choice - Index of the example selected by the user
%     Plots  - Option to display plots: 1.- Yes   2.- No
%     ok     - Configuration acceptance indicator:
%              true.- Accept    false.- Cancel


	%                         Variable initialization                       %
	% ─────────────────────────────────────────────────────────────────────────
	choice = [];
    Plots = 1;
	ok = false;
	% ─────────────────────────────────────────────────────────────────────────


	%                           Window dimensions                           %
	% ─────────────────────────────────────────────────────────────────────────
	w = 420;
	h = 200;
	% ─────────────────────────────────────────────────────────────────────────


	%                            Window creation                             %
	% ─────────────────────────────────────────────────────────────────────────
	fig = dialog( 'Name', 'Initial configuration', ...
		'Position', [0 0 w h], 'WindowStyle', 'modal');

	% Center the window
	movegui(fig, 'center');
	% ─────────────────────────────────────────────────────────────────────────


	%                          Example selection                            %
	% ─────────────────────────────────────────────────────────────────────────
	% Example-list title
	uicontrol(fig, 'Style', 'text', 'String', 'Choose an example:', ...
		'HorizontalAlignment', 'left', 'Position', [30 160 200 25]);

	% List of available examples
	example_list = uicontrol(fig, 'Style', 'listbox', ...
		'String', options, 'Value', 1, 'Position', [30 100 360 66]);
	% ─────────────────────────────────────────────────────────────────────────


	%                           Display options                             %
	% ─────────────────────────────────────────────────────────────────────────
	% Text for plot selection
	uicontrol(fig, 'Style', 'text', 'String', 'Show plots:', ...
		'HorizontalAlignment', 'left', 'Position', [30 60 110 22]);

	% Selector to show or hide plots
	plots_popup = uicontrol(fig, 'Style', 'popupmenu', 'String', {'1.- Yes', '2.- No'}, ...
		'Value', 1, 'Position', [165 60 120 22]);
	% ─────────────────────────────────────────────────────────────────────────


	%                                 Buttons                                %
	% ─────────────────────────────────────────────────────────────────────────
	% Button to accept the configuration
	uicontrol(fig, 'Style', 'pushbutton', 'String', 'Accept', ...
		'Position', [130 20 80 30], 'Callback', @accept_selection);

	% Button to cancel the configuration
	uicontrol(fig, 'Style', 'pushbutton', 'String', 'Cancel', ...
		'Position', [230 20 80 30], 'Callback', @cancel_selection);
	% ─────────────────────────────────────────────────────────────────────────


	%                       Wait and data retrieval                         %
	% ─────────────────────────────────────────────────────────────────────────
	% Wait for the user's response
	uiwait(fig);

	% Close the window if it is still open
	if isvalid(fig)
		delete(fig);
	end
	% ─────────────────────────────────────────────────────────────────────────


	%                         Function to accept                            %
	% ─────────────────────────────────────────────────────────────────────────
	function accept_selection(~, ~)
		choice = get(example_list, 'Value');
		Plots = get(plots_popup, 'Value');
		ok = true;
		uiresume(fig);
	end
	% ─────────────────────────────────────────────────────────────────────────


	%                         Function to cancel                            %
	% ─────────────────────────────────────────────────────────────────────────
	function cancel_selection(~, ~)
		choice = [];
		Plots = [];
		ok = false;
		uiresume(fig);
	end
	% ─────────────────────────────────────────────────────────────────────────

end