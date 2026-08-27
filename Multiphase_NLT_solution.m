function [V_NLT] = Multiphase_NLT_solution(Data)
%{
Function that solves the frequency-dependent cascaded Pi circuit using the Numerical Laplace Transform

Inputs:
    Data - Structure with the following parameters:
    Data(1)   = source_type  % Source_type
    Data(2)   = dt           % Time step (s)
    Data(3)   = t_final      % Observation time (s)
    Data(4)   = N            % Number of samples
    Data(5)   = R_in         % Input resistance (Ohm)
    Data(6)   = L_in         % Input inductance (H)
    Data(7)   = R_L          % Load resistance (Ohm)
    Ll        = l            % Transmission-line length (m)

Outputs:
    V3th            - Voltage at the receiving end of the transmission line
    t_nlt           - Time vector 
    simulation_time - Total solution time (s)
%}


	%                       Start of the simulation                         %
	% ─────────────────────────────────────────────────────────────────────────
	% Start the execution-time counter
	tic;
	% ─────────────────────────────────────────────────────────────────────────
	
	%                          Data extraction                            %
	% ─────────────────────────────────────────────────────────────────────────
	source_type  = Data.Source_type;  % Source_type
	dt      = Data.dt;      % Time step
	T       = Data.t_final; % Final time
	N       = Data.N;       % Number of samples
	R_in    = Data.R_in;    % Input resistance
	L_in    = Data.L_in;    % Input inductance
	R_L     = Data.R_L;     % Load resistance
	Angles  = Data.Angles;  % Phase-angle shifts
	General_angle = Data.General_angle; % General phase-angle shift
	Vac     = Data.Vac;     % Amplitude of source_type
	Vdc     = Data.Vdc;     % Amplitude of source_type DC

	% Numerical NLT options
	if isfield(Data,'Window')
		window = Data.Window;
	else
		window = 3; % Blackman
	end

	if isfield(Data,'Damping')
		damping = Data.Damping;
	else
		damping = 1; % c = 2*dw
	end
	% ─────────────────────────────────────────────────────────────────────────

	%                        NLT discretization                           %
	% ─────────────────────────────────────────────────────────────────────────
	dw = pi/T; % Frequency step, w

	% Factor de damping
	if damping == 1
		c = 2*dw;        % Damping
	elseif damping == 2
		c = log(N^2)/T;  % Damping
	else
		error('Invalid NLT damping option.');
	end

	m = 1:2:2*N;                             % Odd samples
	s = c + 1i*m*dw;                         % Complex frequency, rad/s
	n = 0:N-1;                               % Samples used to calculate Cn (exp(c))
	Cn = (N*2*dw/pi)*exp(c*dt+1i*pi/N).^n;   % Calculation of Cn

	% Window used for the numerical inversion
    if window == 1
		sigma = 0.5*(1+cos(0.5*pi*m/N));                     % Hanning
	elseif window == 2
		sigma = log(2*dw);                                   % Wilcox
	elseif window == 3
		sigma = 0.42 + 0.5*cos(0.5*pi*m/N)+0.08*cos(pi*m/N); % Blackman
	elseif window == 4
		sigma = 1-abs(m/N).^2;                               % Riesz
	elseif window == 5
		sigma = sin(0.5*pi*m/N)./(0.5*pi*m/N);               % Lanczos
	else
		error('Invalid NLT window option.');
    end
	% ─────────────────────────────────────────────────────────────────────────


	%                          NLT parameters                            %
	% ─────────────────────────────────────────────────────────────────────────
	NLT_fitting_data.s  = s;
	NLT_fitting_data.w  = imag(s);
	NLT_fitting_data.f  = NLT_fitting_data.w/(2*pi);
	NLT_fitting_data.Ns = N;
	Aux.Plots = 2;
	load('Parameter_calculation.mat','MTL_data');
	[~, Zs, Ys, Lo] = Three_phase_parameter_calculation(NLT_fitting_data, MTL_data, Aux);
	[~, Nf, N] = size(Zs);
	% ─────────────────────────────────────────────────────────────────────────
	

	%                       Clarke transformation                        %
	% ─────────────────────────────────────────────────────────────────────────
	if Nf == 3

		% Clarke transformation matrix for a three-phase transmission line
		Tc  = (1/sqrt(3))*[   sqrt(2)          0          1;
							-1/sqrt(2)  sqrt(3)/sqrt(2)   1;
							-1/sqrt(2) -sqrt(3)/sqrt(2)   1];

		% Inverse Clarke transformation matrix (T^-1 = T.')
		Tci = Tc.';

	elseif Nf == 6

		% Clarke transformation matrix for a double-circuit transmission line
		Tc = [-1/sqrt(6)  1/sqrt(2)   1/sqrt(6)   1/sqrt(6)        0            0;
			   2/sqrt(6)      0       1/sqrt(6)   1/sqrt(6)        0            0;
			  -1/sqrt(6) -1/sqrt(2)   1/sqrt(6)   1/sqrt(6)        0            0;
				   0          0      -1/sqrt(6)   1/sqrt(6)    1/sqrt(2)    -1/sqrt(6);
				   0          0      -1/sqrt(6)   1/sqrt(6)        0         2/sqrt(6);
				   0          0      -1/sqrt(6)   1/sqrt(6)   -1/sqrt(2)    -1/sqrt(6)];

		% Inverse Clarke transformation matrix for a double-circuit line
		Tci = [-1/sqrt(6)  2/sqrt(6)  -1/sqrt(6)       0            0            0;
				1/sqrt(2)      0      -1/sqrt(2)       0            0            0;
				1/sqrt(6)  1/sqrt(6)   1/sqrt(6)  -1/sqrt(6)   -1/sqrt(6)   -1/sqrt(6);
				1/sqrt(6)  1/sqrt(6)   1/sqrt(6)   1/sqrt(6)    1/sqrt(6)    1/sqrt(6);
				   0           0           0       1/sqrt(2)        0       -1/sqrt(2);
				   0           0           0      -1/sqrt(6)    2/sqrt(6)   -1/sqrt(6)];

	else
		error('The NLT is defined only for 3x3 or 6x6 parameter matrices.');
	end
	% ─────────────────────────────────────────────────────────────────────────


	%               Line model in the frequency domain                 %
	% ─────────────────────────────────────────────────────────────────────────
	for n = 1:N

		ZY = Zs(:,:,n)*Ys(:,:,n);         % Product Z*Y.
	
		% General alternative based on eigenvectors and eigenvalues.
		%[M,L] = eig(A); 
	
		% Using Clarke for all samples. 
		L = Tci*ZY*Tc;
		M = Tc;

		% Calculation of Gamma, coth, and csch using identities
		Gamma = M*sqrt(L)*M^(-1);    % Propagation constant
		Yc = Zs(:,:,n)^(-1)*Gamma;   % Characteristic admittance
		Yn(:,:,n) = Yc;
		ide1 = (eye(Nf) + exp(-2*sqrt(L)*Lo))./(eye(Nf) - exp(-2*sqrt(L)*Lo));
		cotanh = M*diag(diag(ide1))*M^(-1);
		% Alternative expression for parameter A(s) using coth.
		%As(:,:,n) = cotanh*Yc(:,:,n);       % Parameter A(s)

		ide2 = 2./(exp(sqrt(L)*Lo) - exp(-sqrt(L)*Lo));
		cosech = M*diag(diag(ide2))*M^(-1); 
		% Alternative expression for parameter B(s) using csch.
		%Bs(:,:,n) = cosech*Yc(:,:,n);       % Parameter B(s)

		% Using the propagation function
		H  = M*diag(diag(exp(-sqrt(L)*Lo)))*M^(-1);
		H2 = M*diag(diag(exp(-2*sqrt(L)*Lo)))*M^(-1);

		As(:,:,n) = ((eye(Nf) - H2)^(-1)*(H2 + eye(Nf)))*Yc;
		Bs(:,:,n) = ((eye(Nf) - H2)^(-1)*(2*H))*Yc;

	end
	% ─────────────────────────────────────────────────────────────────────────


	%                                 Source_type                             %
	% ─────────────────────────────────────────────────────────────────────────
	fo = 60;          % Frequency, Hz
	wo = 2*pi*fo;     % Angular frequency, rad/s
	Fc = 0;           % Time-shift factor
	Is = zeros(Nf,N); % Equivalent source current

	if source_type == 1 % Source_type AC

		angle = (Angles + General_angle) * (pi/180) - pi/2;

		for k = 1:Nf
			Is(k,:) = Vac(k) * ( (exp(-Fc*s)) .* (s*cos(angle(k)) - wo*sin(angle(k))) ./ (wo.^2 + s.^2) ) / R_in;
		end

	elseif source_type == 2 % Source_type DC

		for k = 1:Nf
			Is(k,:) = Vdc(k) * exp(-Fc*s) ./ (s * R_in);
		end

	else
		error('Invalid source type. Use 1 for AC or 2 for DC.');
	end
	% ─────────────────────────────────────────────────────────────────────────


	%                        Additional variables                        %
	% ─────────────────────────────────────────────────────────────────────────
	% Solution using the transfer function.

	% V1 = zeros(Nf, N);
	% V2 = zeros(Nf, N); 
	% V3 = zeros(Nf, N); 
	% I1 = zeros(Nf, N);
	% I2 = zeros(Nf, N);
	% I3 = zeros(Nf, N); 
	% Ir = zeros(Nf, N);
	% If = zeros(Nf, N);
	% ─────────────────────────────────────────────────────────────────────────

	%                  Frequency-domain system solution                 %
	% ─────────────────────────────────────────────────────────────────────────
	V3 = zeros(Nf,N);

	for k = 1:N
	
		A = (eye(Nf)*1/(L_in))*1/s(k) + eye(Nf)*(1/R_in);
		B = As(:,:,k) + (eye(Nf)*1/(L_in))*1/s(k);
		C = As(:,:,k) + eye(Nf)*1/R_L;
		D = -(eye(Nf)*1/(L_in))*1/s(k);
		E = -Bs(:,:,k);
		Deti = (C*D^2 + A*E^2 - A*B*C)^(-1);

		% Voltages
		% V1(:,k) = (E^2-B*C)*(Deti)*Is(:,k);
		% V2(:,k) = (C*D)*(Deti)*Is(:,k);
		V3(:,k) = (-D*E)*(Deti)*Is(:,k);

		% Currents associated with the nodes and branches
		% I1(:,k) = (eye(Nf)*1/(L_in))*1/s(k)*(V1(:,k) - V2(:,k));
		% I2(:,k) =   As(:,:,k)*V2(:,k) - Bs(:,:,k)*V3(:,k);
		% I3(:,k) = - As(:,:,k)*V3(:,k) + Bs(:,:,k)*V2(:,k);
		% Ir(:,k) = V1(:,k)*(1/R_in);       
		% If(:,k) = V2(:,k)*(1/R_L);  

	end
	% ─────────────────────────────────────────────────────────────────────────


	%                 Time-domain transformation using NLT              %
	% ─────────────────────────────────────────────────────────────────────────

	% The following expressions allow the internal voltages and currents.
    % Currently, only V3 is obtained.

	% V1a(1,:) = real(Cn.*ifft(V1(1,:).*sigma));  
	% V1b(1,:) = real(Cn.*ifft(V1(2,:).*sigma));     
	% V1c(1,:) = real(Cn.*ifft(V1(3,:).*sigma)); 

	% V2a(1,:) = real(Cn.*ifft(V2(1,:).*sigma));  
	% V2b(1,:) = real(Cn.*ifft(V2(2,:).*sigma));     
	% V2c(1,:) = real(Cn.*ifft(V2(3,:).*sigma)); 

	% V3a(1,:) = real(Cn.*ifft(V3(1,:).*sigma));  
	% V3b(1,:) = real(Cn.*ifft(V3(2,:).*sigma));     
	% V3c(1,:) = real(Cn.*ifft(V3(3,:).*sigma)); 

	% I1a(1,:) = real(Cn.*ifft(I1(1,:).*sigma));  
	% I1b(1,:) = real(Cn.*ifft(I1(2,:).*sigma));     
	% I1c(1,:) = real(Cn.*ifft(I1(3,:).*sigma)); 

	% I2a(1,:) = real(Cn.*ifft(I2(1,:).*sigma));  
	% I2b(1,:) = real(Cn.*ifft(I2(2,:).*sigma));     
	% I2c(1,:) = real(Cn.*ifft(I2(3,:).*sigma)); 

	% I3a(1,:) = real(Cn.*ifft(I3(1,:).*sigma));  
	% I3b(1,:) = real(Cn.*ifft(I3(2,:).*sigma));     
	% I3c(1,:) = real(Cn.*ifft(I3(3,:).*sigma));
	% ─────────────────────────────────────────────────────────────────────────

	%                   Output voltage and execution time               %
	% ─────────────────────────────────────────────────────────────────────────
	simulation_time = toc;
	V_NLT = zeros(Nf,N);

	for k = 1:Nf
		V_NLT(k,:) = real(Cn.*ifft(V3(k,:).*sigma));
	end

	fprintf('NLT simulation completed in:               %.4f seconds\n', simulation_time);
	% ─────────────────────────────────────────────────────────────────────────
end