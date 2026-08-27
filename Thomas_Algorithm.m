function x = Thomas_Algorithm(A, d)

% This function solves a linear system whose coefficient matrix
% is tridiagonal using the Thomas algorithm (TDMA -
% Tridiagonal Matrix Algorithm).

% Inputs:
%     A - Tridiagonal coefficient matrix of size n x n
%     d - Right-hand-side vector of size n x 1

% Outputs:
%     x - Solution vector of the system A*x = d


	n = length(d);
	x = zeros(n, 1);


	%                       Extraction of A diagonals                       %
	% ─────────────────────────────────────────────────────────────────────────
	a = diag(A, -1); % Lower
	b = diag(A, 0);  % Main
	c = diag(A, 1);  % Upper
	% ─────────────────────────────────────────────────────────────────────────


	%                   Forward-elimination initialization                  %
	% ─────────────────────────────────────────────────────────────────────────
	c(1) = c(1) / b(1);
	d(1) = d(1) / b(1);
	% ─────────────────────────────────────────────────────────────────────────


	%                          Forward elimination                          %
	% ─────────────────────────────────────────────────────────────────────────
	for i = 2:n-1
		denominator = b(i) - a(i-1) * c(i-1);
		c(i) = c(i) / denominator;
		d(i) = (d(i) - a(i-1) * d(i-1)) / denominator;
	end
	d(n) = (d(n) - a(n-1) * d(n-1)) / (b(n) - a(n-1) * c(n-1));
	% ─────────────────────────────────────────────────────────────────────────


	%                           Back substitution                           %
	% ─────────────────────────────────────────────────────────────────────────
	x(n) = d(n);
	for i = n-1:-1:1
		x(i) = d(i) - c(i) * x(i+1);
	end
	% ─────────────────────────────────────────────────────────────────────────

end