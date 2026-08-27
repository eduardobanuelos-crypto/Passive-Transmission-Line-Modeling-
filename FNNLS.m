function [x]=FNNLS(A,b)  % X = inv(A'*A)*A'*b

% This function solves the Non-Negative Least Squares problem
% using the Fast Non-Negative Least Squares (FNNLS) algorithm,
% developed from the Lawson and Hanson method. The procedure
% precomputes the products A^T A and A^T b to reduce the computational
% cost, and uses active and passive sets to impose the
% non-negativity constraint on the solution variables.

% The code is obtained from: J. S. Sathyakumar, Analysis of Nonnegative 
% Least Squares Algorithms, Major Subject: Ocean Engineering, Texas A&M 
% University, College Station, TX, USA, May 2021. Page 38.

% Inputs:
%     A - System matrix of the least-squares problem
%     b - Right-hand-side vector

% Outputs:
%     x - Solution of the least-squares problem subject to x >= 0


	%                      Matrix-product precomputation                     %
	% ─────────────────────────────────────────────────────────────────────────
	AtA = A.' * A;
	Atb = A.' * b;
	% ─────────────────────────────────────────────────────────────────────────


	%                         Algorithm initialization                       %
	% ─────────────────────────────────────────────────────────────────────────
	n = size(AtA,2);     % Number of variables
	P = [];              % Passive set: variables currently allowed to be > 0
	N = 1:n;             % Active set: variables currently fixed at 0
	x = zeros(n,1);      % Initial solution
	w = (Atb) - (AtA)*x; % Reduced gradient
	wN = w(N);           % Reduced gradient in the active set
	iters = 0;           % Iteration counter
	% ─────────────────────────────────────────────────────────────────────────


	%                            Main iterations                             %
	% ─────────────────────────────────────────────────────────────────────────
	% The process continues while variables remain in the active set and
	% at least one of them has a positive reduced gradient
	tol = 1e-15;
	while(~(isempty(N) || all(wN<=tol)))

		% Active-set variable with the largest reduced gradient
		wmax = max(wN);
		t = find(w==wmax); % Index of the variable entering P
		P = [P t];         % Move variable to the passive set
		N(N==t) = [];      % Remove variable from the active set

		% Least-squares solution restricted to the passive set
		z = zeros(n,1);       % Trial solution
		An = AtA(P,P);
		bn = Atb(P);
		zp = An\bn;           % LS solution in the passive set
		z(N) = zeros(length(N),1);
		z(P) = zp;

		% Inner loop to enforce the non-negativity constraint
		while(min(zp)<=tol)

			% Current solution values in the passive set
			xp = x(P);

			% Ratios used to determine the first variable that reaches zero
			ratio = zeros(length(xp),1);
			for i=1:length(xp)
				ratio(i) = xp(i)/(xp(i)-zp(i));
			end

			qind = zp<=0;                  % Indices that become <= 0
			[alpha,pq] = min(ratio(qind)); % Minimum step that prevents negativity
			q = P(pq);
			x = x + alpha*(z-x);           % Solution update
			temp = P(x(P)==0);
			N = [N temp];                  % Move variables to the active set
			P(P==temp) = [];               % Remove variables from the passive set
			An = AtA(P,P);
			bn = Atb(P);
			zp = An\bn;                    % Recalculate LS in the passive set
			z(N) = zeros(length(N),1);
			z(P) = zp;
		end

		% Update of the solution and reduced gradient
		x=z;
		w = (Atb) - (AtA)*x; % Update the reduced gradient w
		wN = w(N);
		iters = iters +1;
	end
	% ─────────────────────────────────────────────────────────────────────────

end