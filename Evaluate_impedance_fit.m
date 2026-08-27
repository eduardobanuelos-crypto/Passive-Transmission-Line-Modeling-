function [Z_fit,eVF] = Evaluate_impedance_fit(C,P,K,s,Ns,Z,Rdc)

% This function evaluates the rational approximation of an impedance
% function expressed in the form used for the synthesis of a Foster-type
% RL network.

% Inputs:
%     C   - Residues of the rational fitting
%     P   - Poles of the rational fitting
%     K   - Constant term of the rational fitting
%     s   - Complex-frequency samples, rad/s
%     Ns  - Number of frequency samples
%     Z   - Original samples of the impedance function
%     Rdc - Direct-current resistance

% Outputs:
%     Z_fit - Samples of the fitted impedance function
%     eVF   - Absolute error between Z and Z_fit


	%                 Evaluation of the rational approximation                %
	% ─────────────────────────────────────────────────────────────────────────
	Z_fit = zeros(1,Ns);    
	for k = 1:length(P)
		Z_fit = Z_fit + (s.*C(k)./(s - P(k)));
	end 
	Z_fit = Z_fit + Rdc + s.*K;
	% ─────────────────────────────────────────────────────────────────────────


	%                       Calculation of the absolute error                  %
	% ─────────────────────────────────────────────────────────────────────────
	eVF = abs(Z - Z_fit );
	% ─────────────────────────────────────────────────────────────────────────

end