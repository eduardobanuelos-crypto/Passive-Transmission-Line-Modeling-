
% PURPOSE : Calculate the rational approximation f(s)= C*(s*I-A)^(-1)*B + D + s*E
%           where A is diagonal and B is a column of ones (weight function).
%
% The function is programmed as follows
%
%   [C,D]=Residue(Fs,s,Pi,Ns)
%
% INPUTS
%   Fs   : function (column vector) to be fitted. 
%   s    : column vector of frequency points [rad/sec] 
%   Pi   : vector of starting poles
%   Ns   : number of samples
%   Ka   : 1 --> order(numerator)=order(denominator)-1 ('Strictly proper')
%          2 --> order(numerator)=order(denominator)   ('Proper')
%          3 --> order(numerator)=order(denominator)+1 ('Improper')
%
% OUTPUTS
%   C : residues
%   D : constant term
%   E : proportional term
%   A : poles in ascending order
%
function [C,D,E,A]=Residues_nnls(Fs,s,Pi,Ns,Ka,num)

	%   For vectors, SORT(X) sorts the elements of X in ascending order.
	%   For matrices, SORT(X) sorts each column of X in ascending order.
	%   For N-D arrays, SORT(X) sorts the along the first non-singleton dimension of X.
	%   When X is complex, the elements are sorted by ABS(X).  Complex
	%   Matches are further sorted by ANGLE(X).

	% Sort poles in ascending order. First real poles and then complex poles
	Np     = length(Pi);   % Length of the vector that contains the poles
	CPX    = imag(Pi)~=0;  % Put 0 for a real pole and 1 for a complex pole
	rp     = 0;            % Initialize the index for real poles
	cp     = 0;            % Initialize the index for complex poles
	RePole = [];           % Initialize the vector of real poles
	CxPole = [];           % Initialize the vector of complex poles

	% Loop to separate real poles and complex poles
	for k = 1:Np
		if CPX(k) == 0     % Real Pole
			rp = rp + 1;
			RePole(rp) = Pi(k);
		elseif CPX(k) == 1 % Complex pole
			cp = cp + 1;
			CxPole(cp) = Pi(k);
		end
	end

	RePole = sort(RePole);       % Sort real poles
	CxPole = sort(CxPole);       % Sort complex poles
	CxPole = (CxPole.')';        % CxPole=CxPole-2*i*imag(CxPole);
	Lambda = [RePole CxPole];    % Concentrate the full set of starting poles
	I      = diag(ones(1,Np));   % Unit diagonal matrix of ones
	A      = [];                 % Poles
	B      = ones(Ns,1);         % Weight factor (always one for each pole) 
	C      = [];                 % Residues
	D      = zeros(1);           % Initialize the variable to store the constant term (is produced if asympflag=2 or 3)
	E      = zeros(1);           % Initialize the variable to store proportional term (is produced if asympflag=3)

	% Identifies which poles are complex and creates a vector with
	%  0 - for a real pole
	%  1 - for the real part of a complex conjugate poles
	%  2 - for the imaginary part of a complex conjugate poles

	cpx = imag(Lambda)~=0;  % This instruction assigns 0 to a real pole and 1 to a complex pole
	dix = zeros(1,Np);      % Initialize the vector used to identify poles
	if cpx(1)~=0            % If the first pole is complex
		dix(1)=1;           % Put 1 in dix(1) for the real part
		dix(2)=2;           % Put 2 in dix(2) for the imaginary part
		k=3;                % Continue dix from the third position
	else
		k=2;                % If the first pole is real, continue dix from the second position
	end

	% Complete the classification of the poles
	for m=k:Np 
		if cpx(m)~=0         % If the pole is complex
			if dix(m-1)==1
				dix(m)=2;    % If the previous position contains the real part, put 2 to identify the imaginary part
			else
				dix(m)=1;    % Put 1 for the real part of a complex pole
			end
		end
	end

	% Compute the output matrices:
	%   A = Poles (Lambda)
	%   C = Residues
	%   D = Constant term
	%   E = Proportional term

	% This routine builds the Dk matrix equal to 1./(s-Lambda), where Lambda contains the
	% calculated poles; this matrix contains the real poles first
	% followed by the complex poles
	Dk=zeros(Ns,Np);                  
	for m=1:Np
		if dix(m)==0        % Real pole
			Dk(:,m) = B./(s-Lambda(m));
		elseif dix(m)==1    % Complex pole, 1st part
			Dk(:,m) = B./(s-Lambda(m)) + B./(s-Lambda(m)');
		elseif dix(m)==2    % Complex pole, 2st part
			Dk(:,m) = 1i.*B./(s-Lambda(m-1)) - 1i.*B./(s-Lambda(m-1)');
		end
	end 

	% Create the workspace for matrix A and vector b
	AA1=Dk;
	AA2=B.*ones(Ns,1); 
	AA3=B.*s;  

	if Ka == 1
		AA = [AA1];          % Strictly proper rational fitting
	elseif Ka == 2
		AA = [AA1 AA2];      % Proper rational fitting
	elseif Ka == 3
		AA = [AA1 AA2 AA3];  % Improper rational fitting
	else
		disp('Ka need to be 1, 2 or 3')
	end

	bb  = B.*Fs.';

	AAre = real(AA);      % Real part of matrix A
	AAim = imag(AA);      % Imaginary part of matrix A
	bbre = real(bb);      % Real part of matrix b
	bbim = imag(bb);      % Imaginary part of matrix b

	AAn = [AAre; AAim];   % Real and imaginary part of A
	bbn = [bbre; bbim];   % Real and imaginary part of b

	% Solving system  X=inv(A'*A)*A'*b  over-determined system (Ax=b ===> x=A\b)                   
	% if one makes  [Q,R]=qr(A,0) where A=Q*R, one will obtain:
	% X = inv((Q*R)'*Q*R)*(Q*R)'*b, which is X= inv(R'*Q'*Q*R)*R'*Q'*b
	% if Q'*Q=I then X = inv(R'*R)*R'*Q'*b, for this reason one makes A = R and b = Q'*b
	[Q,R]=qr(AAn,0);
	AAn = R;
	bbn = Q.'*bbn;

	[Xmax Ymax] = size(AAn);
	for col=1:Ymax
		Euclidean_norm(col)=norm(AAn(:,col),2);  % Euclidean norm: NORM(V,P) = sum(abs(V).^P)^(1/P).
		AAn(:,col)=AAn(:,col)./Euclidean_norm(col);
	end
 
	%% Solving system  X=inv(A'*A)*A'*b  over-determined system (Ax=b   ===>    x=A\b)  

	if num.NNLS == 1
		Xxn = lsqnonneg(AAn,bbn);
	elseif num.NNLS == 2
		Xxn = FNNLS(AAn,bbn);
	elseif num.NNLS == 3
		Xxn = TNT_NN_method(AAn,bbn);
	elseif num.NNLS == 4
		Xxn = Lawson_Hanson_Algorithm(AAn,bbn);
	elseif num.NNLS == 5
		Xxn = nnls(AAn,bbn);
	end

	X=Xxn./Euclidean_norm.';

	% Put the residues into matrix C
	C=X(1:Np);

	% Make C complex when the residues are complex
	for m=1:Np
		if dix(m)==1
			alpha   = C(m);           % Real part of a complex pole
			beta  = C(m+1);           % Imag part of a complex pole
			C(m)   = alpha + 1i*beta; % The complex pole
			C(m+1) = alpha - 1i*beta; % The conjugate of the previous complex pole
		end
	end

	% Outputs
	if Ka == 1
		A  = Lambda.';   % Poles
		C  = C;          % Residues
		D  = 0;          % Constant term
		E  = 0;          % Proportional term
	elseif Ka == 2
		A  = Lambda.';   % Poles
		C  = C;          % Residues
		D  = X(Np+1);    % Constant term
		E  = 0;          % Proportional term
	elseif Ka == 3
		A  = Lambda.';   % Poles
		C  = C;          % Residues
		D  = X(Np+1);    % Constant term
		E  = X(Np+2);    % Proportional term
	end


end
