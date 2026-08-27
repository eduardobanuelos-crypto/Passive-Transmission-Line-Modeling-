function [Fs_fit,eVF] = Evaluate_rational_fit(C,P,K,E,s,Ns,Fs)

% This function evaluates the rational approximation obtained through
% Vector Fitting from its poles, residues, and constant and proportional
% terms.

% Inputs:
%     C  - Residues of the rational fitting
%     P  - Poles of the rational fitting
%     K  - Constant term of the rational fitting
%     E  - Term proportional to s of the rational fitting
%     s  - Complex-frequency samples
%     Ns - Number of frequency samples
%     Fs - Original samples of the function to be approximated

% Outputs:
%     Fs_fit - Samples of the fitted rational function
%     eVF    - Absolute error between Fs and Fs_fit


	%                 Evaluation of the rational approximation                %
	% ─────────────────────────────────────────────────────────────────────────
	Fs_fit = zeros(1,Ns);    
	for k = 1:length(P)
		Fs_fit = Fs_fit + (C(k)./(s - P(k)));
	end 
	Fs_fit = Fs_fit + K + E.*s;
	% ─────────────────────────────────────────────────────────────────────────

	%                       Calculation of the absolute error                  %
	% ─────────────────────────────────────────────────────────────────────────
	eVF = abs(Fs - Fs_fit );
	% ─────────────────────────────────────────────────────────────────────────

end