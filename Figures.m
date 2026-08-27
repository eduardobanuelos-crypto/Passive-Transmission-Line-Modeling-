function []=Figures()

% This function sets the default graphical configuration
% used in the program figures.

	%                         General figure format                         %
	% ─────────────────────────────────────────────────────────────────────────
	set(groot, 'DefaultFigureWindowStyle','normal');
	set(groot, 'DefaultAxesFontName',     'Times New Roman');
	set(groot, 'DefaultLegendFontName',   'Times New Roman');
	set(groot, 'DefaultFigureColor',      'white');
	set(groot, 'DefaultAxesFontWeight',   'normal');
	set(groot, 'DefaultLegendLocation',   'northeast');
	set(groot, 'DefaultTextInterpreter',  'latex');
	set(groot, 'DefaultAxesXGrid',        'on');
	set(groot, 'DefaultAxesYGrid',        'on');
	set(groot, 'DefaultAxesFontSize',      24);
	set(groot, 'DefaultLineLineWidth',     2.5);
	set(groot, 'DefaultFigureUnits',      'pixels');
	set(groot, 'DefaultFigurePosition',   [0 100 1600 500]);  % width x height (pixels)
	%set(groot, 'DefaultFigureUnits',     'centimeters');
	% set(groot, 'DefaultFigurePosition', [2 2 40 8]);        % width x height (cm)
	set(groot, 'DefaultFigureColormap',    winter);
	% ─────────────────────────────────────────────────────────────────────────

end