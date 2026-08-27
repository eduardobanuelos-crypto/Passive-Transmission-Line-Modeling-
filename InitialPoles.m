function [Ps]=InitialPoles(f,Npol)

% This function generates the set of initial poles used by Vector Fitting.

% Inputs:
%     f    - Frequency vector used to distribute the initial
%            poles
%     Npol - Total number of initial poles
%
% Outputs:
%     Ps - Column vector containing the initial poles


	%                    Determination of the number of poles                %
	% ─────────────────────────────────────────────────────────────────────────
	n_complex_pairs  = fix(Npol/2);           % Number of complex-conjugate initial pole pairs
	odd_remainder = Npol/2 - n_complex_pairs; % Auxiliary variable used to check whether Npol is odd
	has_real_pole  = odd_remainder ~= 0;      % Put  0 - n_complex_pairs initial poles  &  1 - odd initial poles
	% ─────────────────────────────────────────────────────────────────────────


	%                          Pole initialization                          %
	% ─────────────────────────────────────────────────────────────────────────
	% Set a real pole when has_real_pole == 1
	if has_real_pole == 0    % Even number of initial poles
		poles = [];
    else                     % Odd number of initial poles
		poles = [(max(f)-min(f))/2];
	end
	% ─────────────────────────────────────────────────────────────────────────


	%                   Generation of complex-conjugate poles               %
	% ─────────────────────────────────────────────────────────────────────────
	% Set the complex-conjugate initial poles
	beta = linspace(min(f),max(f),n_complex_pairs);
	for n=1:length(beta)
	  alpha=-beta(n)*1e-2;
	  poles=[poles (alpha-1j*beta(n)) (alpha+1j*beta(n)) ]; 
	end
	% ─────────────────────────────────────────────────────────────────────────


	%                          Initial-pole vector                          %
	% ─────────────────────────────────────────────────────────────────────────
	Ps = poles.';  % Column vector of the initial poles
	% ─────────────────────────────────────────────────────────────────────────

end