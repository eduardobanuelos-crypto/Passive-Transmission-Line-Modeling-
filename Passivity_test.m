function Passivity_test()

% This function evaluates the passivity of the transmission-line model
% represented in state space. For each available mode,
% the transfer matrix is obtained in the frequency domain and
% its Hermitian part is calculated. The passivity condition is verified
% through the eigenvalues of this matrix, which must remain non-negative
% over the analyzed frequency range.


	%                            Matrix loading                             %
	% ─────────────────────────────────────────────────────────────────────────
	load("Passivity_matrices.mat",'Passivity_matrices','Data');
	Matrices = Passivity_matrices;

	f = Data.Fitting_data.f;
	s = Data.Fitting_data.s;
	N = length(s);

	% Number of modes available in the modal representation
	Nf = numel(Matrices.SS_matrix);

	if Nf ~= 3 && Nf ~= 6
		error('The passivity test is implemented only for 3x3 or 6x6 cases.');
	end
	% ─────────────────────────────────────────────────────────────────────────


	%                         Alpha and beta modes                          %
	% ─────────────────────────────────────────────────────────────────────────
	A_alpha = Matrices.SS_matrix(1).A;
	B_alpha = Matrices.SS_m(1).B;
	C_alpha = Matrices.SS_m2(1).C;
	D_alpha = Matrices.SS_m3(1).D;

	[H, ZH, lam] = calculate_realization(A_alpha, B_alpha, C_alpha, D_alpha, s, N);
	% ─────────────────────────────────────────────────────────────────────────


	%                            First zero mode                            %
	% ─────────────────────────────────────────────────────────────────────────
	A_zero = Matrices.SS_matrix(3).A;
	B_zero = Matrices.SS_m(3).B;
	C_zero = Matrices.SS_m2(3).C;
	D_zero = Matrices.SS_m3(3).D;

	[H_zero, ZH_zero, lam_zero] = calculate_realization(A_zero, B_zero, C_zero, D_zero, s, N);
	% ─────────────────────────────────────────────────────────────────────────


	%                       Second zero mode: 6x6 case                     %
	% ─────────────────────────────────────────────────────────────────────────
	if Nf == 6

		A_zero_2 = Matrices.SS_matrix(4).A;
		B_zero_2 = Matrices.SS_m(4).B;
		C_zero_2 = Matrices.SS_m2(4).C;
		D_zero_2 = Matrices.SS_m3(4).D;

		[H_zero_2, ZH_zero_2, lam_zero_2] = ...
			calculate_realization(A_zero_2, B_zero_2, C_zero_2, D_zero_2, s, N);

	end
	% ─────────────────────────────────────────────────────────────────────────


	%                           Plot of H elements                          %
	% ─────────────────────────────────────────────────────────────────────────
	figure;
	map = 'turbo';
	colormap(map);

	for n = 1:3
		for m = 1:3

			X = f;

			% Transfer matrix of the alpha and beta modes
			Y = squeeze(abs(H(n,m,:)))';
			C = log10(max(Y,realmin));

			surface([X;X], [Y;Y], [zeros(size(X));zeros(size(X))], ...
					[C;C], 'FaceColor','none', 'EdgeColor','interp', ...
					'HandleVisibility','off', 'LineWidth',2.5);
			hold on;

			% Transfer matrix of the first zero mode
			Y2 = squeeze(abs(H_zero(n,m,:)))';
			C2 = log10(max(Y2,realmin));

			surface([X;X], [Y2;Y2], [zeros(size(X));zeros(size(X))], ...
					[C2;C2], 'FaceColor','none', 'EdgeColor','interp', ...
					'HandleVisibility','off', 'LineWidth',2.5);

			% Transfer matrix of the second zero mode for the 6x6 case
			if Nf == 6
				Y3 = squeeze(abs(H_zero_2(n,m,:)))';
				C3 = log10(max(Y3,realmin));

				surface([X;X], [Y3;Y3], [zeros(size(X));zeros(size(X))], ...
						[C3;C3], 'FaceColor','none', 'EdgeColor','interp', ...
						'HandleVisibility','off', 'LineWidth',2.5);
			end

		end
	end

	set(gca,'XScale','log');
	set(gca,'YScale','log');
	grid on;
	xlabel('Frequency [Hz]');
	ylabel('Magnitude');

	colorbar;
	clim([-2 4]);

	xlim([min(f) max(f)]);
	ylim([1e-11 1e5]);
	yticks([1e-11 1e-7 1e-3 1e1 1e5]);
	set(gcf,'pos',[150 300 1600 500]);
	% ─────────────────────────────────────────────────────────────────────────


	%                           Eigenvalue plot                             %
	% ─────────────────────────────────────────────────────────────────────────
	figure;
	map = 'turbo';
	colormap(map);

	for n = 1:3

		X = f;

		% Eigenvalues of the alpha and beta modes
		Y = squeeze(lam(n,1,:))';
		C = Y;

		surface([X;X], [Y;Y], [zeros(size(X));zeros(size(X))], ...
				[C;C], 'FaceColor','none', 'EdgeColor','interp', ...
				'HandleVisibility','off', 'LineWidth',2.5);
		hold on;

		% Eigenvalues of the first zero mode
		Y2 = squeeze(lam_zero(n,1,:))';
		C2 = Y2;

		surface([X;X], [Y2;Y2], [zeros(size(X));zeros(size(X))], ...
				[C2;C2], 'FaceColor','none', 'EdgeColor','interp', ...
				'HandleVisibility','off', 'LineWidth',2.5);

		% Eigenvalues of the second zero mode for the 6x6 case
		if Nf == 6
			Y3 = squeeze(lam_zero_2(n,1,:))';
			C3 = Y3;

			surface([X;X], [Y3;Y3], [zeros(size(X));zeros(size(X))], ...
					[C3;C3], 'FaceColor','none', 'EdgeColor','interp', ...
					'HandleVisibility','off', 'LineWidth',2.5);
		end

	end

	set(gca,'XScale','log');
	grid on;
	xlabel('Frequency [Hz]');
	ylabel('Amplitude');

	colorbar;
	clim([-1 5]);

	xlim([min(f) max(f)]);
	ylim([-1 5]);
	yticks([-1 0 1 2 3 4 5]);
	set(gcf,'pos',[150 300 1600 500]);
	% ─────────────────────────────────────────────────────────────────────────

end


function [H, ZH, lam] = calculate_realization(A, B, C, D, s, N)

	n = size(A,1);
	I_n = eye(n);

	H   = zeros(3,3,N);
	ZH  = zeros(3,3,N);
	lam = zeros(3,1,N);

    for k = 1:N
		H(:,:,k)   = C * ((s(1,k)*I_n - A) \ B) + D;
		ZH(:,:,k)  = (H(:,:,k) + H(:,:,k)')/2;
		lam(:,:,k) = eig(ZH(:,:,k));
    end

end