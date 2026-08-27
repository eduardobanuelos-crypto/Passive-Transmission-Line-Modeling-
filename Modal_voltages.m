function [Modal_voltages] = Modal_voltages(source_type, t, Voltage_values)

% This function generates the source voltages in the phase domain
% and transforms them to the modal domain using the Clarke
% transform.

% Inputs:
%     source_type    - Source type: 1.- AC    2.- DC
%     t              - Simulation time vector, s
%     Voltage_values - Structure with the magnitudes and angles of the
%                      source voltages

% Outputs:
%     Modal_voltages - Structure with the source voltages in the
%                      modal domain for the considered configuration


	%                            Data extraction                            %
	% ─────────────────────────────────────────────────────────────────────────
	Vac           = Voltage_values.Vac;           % AC-source voltages
	Angles        = Voltage_values.Angles;        % Phase-angle shifts
	General_angle = Voltage_values.General_angle; % General phase-angle shift
	Vdc           = Voltage_values.Vdc;           % DC-source voltages
	Nf            = length(Vac);                  % Number of source phases/conductors
	% ─────────────────────────────────────────────────────────────────────────


	%                            Source voltages                            %
	% ─────────────────────────────────────────────────────────────────────────
	if source_type == 1 % AC source

		V_source = zeros(Nf,length(t));

		for k = 1:Nf
			V_source(k,:) = Vac(k,1) * sin(2*pi*60*t + Angles(k,1)*(pi/180) + General_angle*(pi/180));
		end

	else % DC source

		V_source = zeros(Nf,length(t));

		for k = 1:Nf
			V_source(k,:) = Vdc(k,1) * ones(1,length(t));
		end

	end
	% ─────────────────────────────────────────────────────────────────────────


	%                       Clarke transformation                         %
	% ─────────────────────────────────────────────────────────────────────────
	if Nf == 3

		% Clarke transform for a three-phase transmission line
		TClarke = sqrt(2/3) * [    1,       -1/2,       -1/2;
								   0,     sqrt(3)/2,  -sqrt(3)/2;
							   1/sqrt(2), 1/sqrt(2),   1/sqrt(2)];

	elseif Nf == 6

		% Clarke transform for a double-circuit three-phase transmission line
		TClarke = [-1/sqrt(6), 2/sqrt(6),  -1/sqrt(6),     0,            0,           0;
					1/sqrt(2),     0,      -1/sqrt(2),      0,            0,           0;
					1/sqrt(6), 1/sqrt(6),   1/sqrt(6),  -1/sqrt(6),  -1/sqrt(6),  -1/sqrt(6);
					1/sqrt(6), 1/sqrt(6),   1/sqrt(6),   1/sqrt(6),   1/sqrt(6),   1/sqrt(6);
					   0,          0,           0,       1/sqrt(2),       0,      -1/sqrt(2);
					   0,          0,           0,      -1/sqrt(6),   2/sqrt(6),  -1/sqrt(6)];

	else
		error('Modal_voltages is implemented only for 3- or 6-phase cases.');
	end
	% ─────────────────────────────────────────────────────────────────────────


	%                      Transformation to modal domain                   %
	% ─────────────────────────────────────────────────────────────────────────
	V_alpha_beta_0 = TClarke * V_source;
	% ─────────────────────────────────────────────────────────────────────────


	%                       Voltages in the modal domain                    %
	% ─────────────────────────────────────────────────────────────────────────
	if Nf == 3

		Modal_voltages.V_alpha = V_alpha_beta_0(1,:);
		Modal_voltages.V_beta  = V_alpha_beta_0(2,:);
		Modal_voltages.V_zero  = V_alpha_beta_0(3,:);

	elseif Nf == 6

		Modal_voltages.V_alpha  = V_alpha_beta_0(1,:);
		Modal_voltages.V_beta   = V_alpha_beta_0(2,:);
		Modal_voltages.V_zero   = V_alpha_beta_0(3,:);
		Modal_voltages.V_zero2  = V_alpha_beta_0(4,:);
		Modal_voltages.V_beta2  = V_alpha_beta_0(5,:);
		Modal_voltages.V_alpha2 = V_alpha_beta_0(6,:);

	end
	% ─────────────────────────────────────────────────────────────────────────

end