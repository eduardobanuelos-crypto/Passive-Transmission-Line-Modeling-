function [x] = Lawson_Hanson_Algorithm(A,b)

% This function solves the Non-Negative Least Squares problem
% using the Lawson and Hanson algorithm. The procedure
% uses active and passive sets to identify the variables that
% remain constrained to zero and those that may take positive
% values, iteratively solving least-squares problems until
% the non-negativity constraint is satisfied.

% The code is obtained from: J. S. Sathyakumar, Analysis of Nonnegative Least Squares Algorithms, Major Sub-
% ject: Ocean Engineering, Texas A&M University, College Station, TX, USA, May 2021. Page 37.

% Inputs:
%     A - System matrix of the least-squares problem
%     b - Right-hand-side vector

% Outputs:
%     x - Solution of the least-squares problem subject to x >= 0


	%                         Algorithm initialization                       %
	% ─────────────────────────────────────────────────────────────────────────
	n = size(A,2);
	P = [];
	N = 1:n;
	x = zeros(n,1);
	w = A'*(b - A*x); % Negative gradient
	wN = w(N);
	iters = 0;
	% ─────────────────────────────────────────────────────────────────────────


	%                            Main iterations                             %
	% ─────────────────────────────────────────────────────────────────────────
	while(~(isempty(N) || all(wN<=1e-12)))

		% Transfer of variables between the active and passive sets
		wmax = max(wN);
		t = find(w==wmax);
		P = [P t];
		N(N==t) = [];

		% Least-squares solution in the passive set
		z = zeros(n,1);
		clear("Ap");
		Ap = A(:,P);
		zp = Ap\b;
		z(N) = zeros(length(N),1);
		z(P) = zp;

		% Inner loop to enforce the non-negativity constraint
		while(min(zp)<=0)
			xp = x(P);
			ratio = zeros(length(xp),1);

			% Calculation of the ratios used to determine the admissible step
			for i=1:length(xp)
				ratio(i) = xp(i)/(xp(i)-zp(i));
			end

			qind = zp<=0;
			[alpha,pq] = min(ratio(qind));
			q = P(pq);
			x = x + alpha*(z-x);
			temp = P(x(P)==0);
			N = [N temp];
			P(P==temp) = [];
			clear("Ap");
			Ap = A(:,P);
			zp = Ap\b;
			z(N) = zeros(length(N),1);
			z(P) = zp;
		end

		% Update of the solution and reduced gradient
		x=z;
		w = A'*(b - A*x);
		wN = w(N);
		iters = iters + 1;
	end
	% ─────────────────────────────────────────────────────────────────────────

end