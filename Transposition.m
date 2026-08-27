function [Z_transposed, Y_transposed] = Transposition(Z_tot, Y_tot, n_conds, N_samples)

% This function obtains the impedance and admittance matrices of
% an ideally transposed transmission line. The implementation
% considers three-phase and double-circuit lines.

% Inputs:
%     Z_tot     - Line impedance matrix in the phase domain
%     Y_tot     - Line admittance matrix in the phase domain
%     n_conds   - Number of phase conductors considered
%     N_samples - Number of frequency samples

% Outputs:
%     Z_transposed - Impedance matrix of the ideally transposed line
%     Y_transposed - Admittance matrix of the ideally transposed line


	%                         Variable initialization                       %
	% ─────────────────────────────────────────────────────────────────────────
	Nf = n_conds;
	Z_transposed = zeros(Nf, Nf, N_samples);
	Y_transposed = zeros(Nf, Nf, N_samples);
	% ─────────────────────────────────────────────────────────────────────────


	if Nf == 3

		%                   Three-phase line transposition                      %
		% ─────────────────────────────────────────────────────────────────────────
		% Auxiliary impedance matrices
		Z1 = zeros(Nf, Nf, N_samples);
		Z2 = zeros(Nf, Nf, N_samples);
		Z3 = zeros(Nf, Nf, N_samples);

		% Auxiliary admittance matrices
		Y1 = zeros(Nf, Nf, N_samples);
		Y2 = zeros(Nf, Nf, N_samples);
		Y3 = zeros(Nf, Nf, N_samples);

		for k = 1:N_samples
			% First impedance-transposition position
			Z1(:,:,k) = Z_tot(:,:,k);

			% Second transposition position: A->B, B->C, C->A
			Z2(:,:,k) = [Z_tot(2,2,k), Z_tot(2,3,k), Z_tot(2,1,k);
						 Z_tot(3,2,k), Z_tot(3,3,k), Z_tot(3,1,k);
						 Z_tot(1,2,k), Z_tot(1,3,k), Z_tot(1,1,k)];

			% Third transposition position: A->C, B->A, C->B
			Z3(:,:,k) = [Z_tot(3,3,k), Z_tot(3,1,k), Z_tot(3,2,k);
						 Z_tot(1,3,k), Z_tot(1,1,k), Z_tot(1,2,k);
						 Z_tot(2,3,k), Z_tot(2,1,k), Z_tot(2,2,k)];

			% Averaged transposed impedance matrix
			Z_transposed(:,:,k) = (1/3) .* (Z1(:,:,k) + Z2(:,:,k) + Z3(:,:,k));

			% First admittance-transposition position
			Y1(:,:,k) = Y_tot(:,:,k);

			% Second transposition position: A->B, B->C, C->A
			Y2(:,:,k) = [Y_tot(2,2,k), Y_tot(2,3,k), Y_tot(2,1,k);
						 Y_tot(3,2,k), Y_tot(3,3,k), Y_tot(3,1,k);
						 Y_tot(1,2,k), Y_tot(1,3,k), Y_tot(1,1,k)];

			% Third transposition position: A->C, B->A, C->B
			Y3(:,:,k) = [Y_tot(3,3,k), Y_tot(3,1,k), Y_tot(3,2,k);
						 Y_tot(1,3,k), Y_tot(1,1,k), Y_tot(1,2,k);
						 Y_tot(2,3,k), Y_tot(2,1,k), Y_tot(2,2,k)];

			% Averaged transposed admittance matrix
			Y_transposed(:,:,k) = (1/3) .* (Y1(:,:,k) + Y2(:,:,k) + Y3(:,:,k));
		end
		% ─────────────────────────────────────────────────────────────────────────

	elseif Nf == 6

		for k = 1:N_samples

			%             Impedance transposition for a double circuit             %
			% ─────────────────────────────────────────────────────────────────────────
			% Average self impedance
			zs = ( Z_tot(1,1,k) + Z_tot(2,2,k) + Z_tot(3,3,k) + ...
				   Z_tot(4,4,k) + Z_tot(5,5,k) + Z_tot(6,6,k) ) / 6;


			% Average coupling within each circuit
			% Circuit I: AI-BI, AI-CI, BI-CI
			% Circuit II: AII-BII, AII-CII, BII-CII

			zm = ( Z_tot(1,2,k) + Z_tot(1,3,k) + Z_tot(2,3,k) + ...
				   Z_tot(4,5,k) + Z_tot(4,6,k) + Z_tot(5,6,k) ) / 6;

			% Average coupling between circuits
			% AI-AII, AI-BII, AI-CII
			% BI-AII, BI-BII, BI-CII
			% CI-AII, CI-BII, CI-CII

			zp = ( Z_tot(1,4,k) + Z_tot(1,5,k) + Z_tot(1,6,k) + ...
				   Z_tot(2,4,k) + Z_tot(2,5,k) + Z_tot(2,6,k) + ...
				   Z_tot(3,4,k) + Z_tot(3,5,k) + Z_tot(3,6,k) ) / 9;


			% Transposed impedance matrix
			Z_transposed(:,:,k) = [ zs  zm  zm  zp  zp  zp;
								    zm  zs  zm  zp  zp  zp;
								    zm  zm  zs  zp  zp  zp;
								    zp  zp  zp  zs  zm  zm;
									zp  zp  zp  zm  zs  zm;
									zp  zp  zp  zm  zm  zs ];
			% ─────────────────────────────────────────────────────────────────────────


			%             Admittance transposition for a double circuit            %
			% ─────────────────────────────────────────────────────────────────────────
			% Average self admittance
			ys = ( Y_tot(1,1,k) + Y_tot(2,2,k) + Y_tot(3,3,k) + ...
				   Y_tot(4,4,k) + Y_tot(5,5,k) + Y_tot(6,6,k) ) / 6;


			% Average coupling within each circuit
			ym = ( Y_tot(1,2,k) + Y_tot(1,3,k) + Y_tot(2,3,k) + ...
				   Y_tot(4,5,k) + Y_tot(4,6,k) + Y_tot(5,6,k) ) / 6;


			% Average coupling between circuits
			yp = ( Y_tot(1,4,k) + Y_tot(1,5,k) + Y_tot(1,6,k) + ...
				   Y_tot(2,4,k) + Y_tot(2,5,k) + Y_tot(2,6,k) + ...
				   Y_tot(3,4,k) + Y_tot(3,5,k) + Y_tot(3,6,k) ) / 9;

			% Transposed admittance matrix
			Y_transposed(:,:,k) = [ ys  ym  ym  yp  yp  yp;
									ym  ys  ym  yp  yp  yp;
									ym  ym  ys  yp  yp  yp;
									yp  yp  yp  ys  ym  ym;
									yp  yp  yp  ym  ys  ym;
									yp  yp  yp  ym  ym  ys ];
		end
		% ─────────────────────────────────────────────────────────────────────────

	end
	
end