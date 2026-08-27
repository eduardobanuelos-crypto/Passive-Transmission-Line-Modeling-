function [Values, Zp_fit] = VF_VFPR(Fitting_data,options,Zp,Ka,R_dc,Ll,Y_shunt,Z, Plot_index)

% This function performs the rational fitting of the series parameters
% of the transmission line using conventional Vector Fitting (VF) and
% Vector Fitting with Real Poles (VFPR). Starting from the initial VF fitting,
% the complex poles are replaced by real poles and the residues are then
% calculated using Non-Negative Least Squares (NNLS).
% Finally, the function obtains the R, L, and C elements of the
% equivalent Pi circuit through Foster synthesis, evaluates the rational
% approximation, and determines the poles and zeros of the fitted function.

% Inputs:
%     Fitting_data - Structure with the parameters used in the fitting
%     options      - Structure with the VFPR, NNLS, and plot options
%     Zp           - Samples of the function used for rational fitting
%     Ka           - Option associated with the asymptotic term used by VF
%     R_dc         - Direct-current resistance per unit length
%     Ll           - Transmission-line length, m
%     Y_shunt      - Shunt admittance used to obtain C
%     Z            - Samples of the original series impedance
%     Plot_index   - Indicator of the mode considered for the plots

% Outputs:
%     Values - Structure with the R, L, C, and G parameters of the Pi circuit
%     Zp_fit - Samples of the rational approximation of Zp
%     Zeros  - Zeros of the rational function obtained with VFPR-NNLS
%     Poles  - Final poles of the rational function after NNLS
%     P_VF   - Poles obtained using conventional Vector Fitting
%     C_VF   - Residues obtained using conventional Vector Fitting


	%                            Data extraction                            %
	% ─────────────────────────────────────────────────────────────────────────
	Ns  = Fitting_data.Ns;
	iterations = Fitting_data.iterations;
	Np  = Fitting_data.Np;
	M   = Fitting_data.M;
	f   = Fitting_data.f;
	w   = Fitting_data.w;
	s   = Fitting_data.s;
	% ─────────────────────────────────────────────────────────────────────────


	%                    Approximation with Vector Fitting                  %
	% ─────────────────────────────────────────────────────────────────────────
	% Initial poles
	[Ps] = InitialPoles(f,Np);

	% Vector Fitting options
	weight=ones(1,Ns); % Vector of weights
	opts.relax=1;      % Use Vector Fitting with a relaxed non-triviality constraint
	opts.stable=1;     % Enforce stable poles
	opts.asymp=Ka;     % Include D and E in the fitting    
	opts.skip_pole=0;  % Do NOT skip pole identification
	opts.skip_res=1;   % Do skip identification of residues (C,D,E) 
	opts.cmplx_ss=1;   % Create a real-only state-space model

	opts.spy1=0;       % No plotting for the first stage of Vector Fitting
	opts.spy2=0;       % Create a magnitude plot for the fitting of f(s) 
	opts.logx=1;       % Use a linear abscissa axis
	opts.logy=1;       % Use a logarithmic ordinate axis 
	opts.errplot=1;    % Include deviation in the magnitude plot
	opts.phaseplot=0;  % Do NOT produce a phase-angle plot
	opts.legend=1;     % Include legends in plots

	% Vector Fitting iterations
	for k = 1:iterations
		if k==iterations 
			opts.skip_res = 0; 
		end
		[SER,Ps,rmserr_VF,fit] = vectfit3(Zp,s,Ps,weight,opts);
	end

	% Parameters obtained using VF
	C_VF = SER.C; % Residues
	D_VF = SER.D; % K
	E_VF = SER.E; % E
	P_VF = Ps;    % Poles
	% ─────────────────────────────────────────────────────────────────────────


	%                            First VFPR stage                            %
	% ─────────────────────────────────────────────────────────────────────────
	C_VFPR = SER.C;
	D_VFPR = SER.D;
	E_VFPR = SER.E;

	% Separation of complex and real poles
	fac = 0.1;
	Real_poles = [];
	Complex_poles = [];
	Pc = [];
	l = 1;
	h = 1;
	cpx  = imag(Ps)~=0;
	for k = 1:1:Np
		if cpx(k) == 0
			Real_poles(h,1) = Ps(k);
			h = h+1;
		end
		if cpx(k) == 1
			Complex_poles(l,1) = Ps(k);
			l = l+1;
		end
	end

	% Replacement of each complex-conjugate pole pair by real poles
    for k = 1:2:length(Complex_poles)
		Wn = abs(Complex_poles(k));
		Pc(k,1) = -(1+fac)*Wn;
		Pc(k+1,1) = -(1-fac)*Wn;
    end

	% Combine the original real poles with the new real poles
	P_VFPR = sort([Real_poles; Pc]);
	% ─────────────────────────────────────────────────────────────────────────


	%                           Second VFPR stage                            %
	% ─────────────────────────────────────────────────────────────────────────
	% Residue calculation using Non-Negative Least Squares
	[C_VFPR,K_VFPR,E_VFPR,P_VFPR] = Residues_nnls(Zp,s.',P_VFPR,Ns,Ka,options);

	% Elimination of poles associated with zero residues
	m = 1; P = [];
	for k = 1:1:length(P_VFPR)   % Remove poles whose residue is zero 
		if C_VFPR(k)~=0
			P(m,1) = P_VFPR(k);
			C(m,1) = C_VFPR(k);
			m = m+1; 
		end
	end

	Eliminated_poles = length(P_VFPR)-length(P);
	% ─────────────────────────────────────────────────────────────────────────


	%                Equivalent Pi-circuit parameter calculation             %
	% ─────────────────────────────────────────────────────────────────────────
	Ro = (R_dc * Ll) / M;                   % Series resistance, Ohm
	Lo = (K_VFPR * Ll) / M;                 % Series inductance, H
	Lo_ATP = 1000 * (K_VFPR * Ll) / M;      % Inductance for ATP, mH

	R_k = (C * Ll) ./ M;                    % Parallel resistances, Ohm
	L_k = ((-C ./ P) * Ll) ./ M;            % Parallel inductances, H
	L_k_ATP = 1000 * ((-C ./ P) * Ll) ./ M; % Inductances for ATP, mH

	Cap = zeros(1,Ns);                      % Capacitance p.u.l., F/m
	for k = 1:1:Ns
		Cap(1,k) = imag(Y_shunt/w);         % Shunt capacitance, F/m
	end
	C_k = (Cap(1) * Ll) / M;                % Total capacitance, F
	C_k_ATP = 1e6 * (Cap(1) * Ll) / M;      % Capacitance for ATP, μF

	G = 0;                                  % Shunt conductance

	% Structure with the R, L, and C parameters of the equivalent circuit
	Values = struct('Ro',Ro,'Rk',R_k,'Lo',Lo,'Lo_ATP',Lo_ATP,'Lk',L_k, ...
		'Lk_ATP',L_k_ATP,'Ck',C_k,'Ck_ATP',C_k_ATP,'G',G);
	% ─────────────────────────────────────────────────────────────────────────


	%                            Fitting evaluation                          %
	% ─────────────────────────────────────────────────────────────────────────
	[Zp_fit,eVF1] = Evaluate_rational_fit(C,P,K_VFPR,E_VFPR,s,Ns,Zp); % Evaluation of the Zp fitting
	[Z_fit, eVF2] = Evaluate_impedance_fit(C,P,K_VFPR,s,Ns,Z,R_dc);   % Evaluation of the Z fitting
	% ─────────────────────────────────────────────────────────────────────────


	%                  Residue-Pole to Zero-Pole conversion                 %
	% ─────────────────────────────────────────────────────────────────────────
	[b,a] = residue(C,P,K_VFPR); % Residue-Pole to Zero-Pole conversion
	Zeros = roots(b);            % Zeros
	Poles = P;                   % Poles
	% ─────────────────────────────────────────────────────────────────────────


	%                                  Plots                                 %
	% ─────────────────────────────────────────────────────────────────────────
	if options.Plots == 1

		% Comparison of the rational fittings
		figure;
		subplot(1,2,1)
		loglog(f,abs(Z),'-','Color',[0 0 0]), hold on
		loglog(f,abs(Z_fit),'--','Color',[0 0 1])
		loglog(f,eVF2,'-','Color',[0 0 1]), hold off
		xlabel({'Frequency [Hz]','(b)'}) 
		ylabel('Magnitude [p.u.l.]')
		legend('Data','VF-RP-NNLS', 'Deviation','Location','southeast'), 

		subplot(1,2,2)
		loglog(f,abs(Zp),'-','Color',[0 0 0]),hold on
		loglog(f,abs(Zp_fit),'--','Color',[0 0 1])
		loglog(f,eVF1,'-','Color',[0 0 1]), hold off
		xlabel({'Frequency [Hz]','(a)'}) 
		ylabel('Magnitude [p.u.l.]')
		legend('Data','VF-RP-NNLS', 'Deviation','Location','southwest')

		% Poles and zeros obtained using VF and VFPR-NNLS
		figure;
		semilogx(real(P_VF./(2*pi)),imag(P_VF./(2*pi)),'x','Color',[0 1 0],'MarkerSize',12), hold on
		semilogx(real(P_VFPR./(2*pi)),imag(P_VFPR./(2*pi)),'x','Color',[0.7 0.7 0.7],'MarkerSize',12);
		semilogx(real(Poles./(2*pi)),imag(Poles./(2*pi)),'x','Color',[0 0 1],'MarkerSize',12);
		semilogx(real(Zeros./(2*pi)),imag(Zeros./(2*pi)),'o','Color',[0.8 0.3 0.2],'MarkerSize',12);
		xlabel({'Real part','(b)'}); ylabel('Imaginary part');
		legend('Poles VF','Eliminated Poles','Poles (VF-RP-NNLS)','Zeros (VF-RP-NNLS)','Orientation','horizontal','Location','northwest');
		if Plot_index == 1
			title(sprintf('N = %d poles with %d eliminated by NNLS in the alpha and beta modes', Np, Eliminated_poles));
		elseif Plot_index == 2
			title(sprintf('N = %d poles with %d eliminated by NNLS in the zero mode', Np, Eliminated_poles));
		elseif Plot_index == 3
			title(sprintf('N = %d poles with %d eliminated by NNLS in the common zero mode', Np, Eliminated_poles));
		end

	end
	% ─────────────────────────────────────────────────────────────────────────

end