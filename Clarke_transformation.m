function [Transformed_impedance, Transformed_admittance] = Clarke_transformation(Z_transposed, Y_transposed, N_samples, n_cnds)

% This function transforms the impedance and admittance matrices
% of a transmission line from the phase domain to the modal domain
% using the Clarke transform.

% Inputs:
%     Z_transposed - Impedance matrix in the phase domain
%     Y_transposed - Admittance matrix in the phase domain
%     N_samples    - Number of frequency samples
%     n_cnds       - Number of phase conductors considered:
%                     3.- Three-phase line    6.- Double circuit

% Outputs:
%     Transformed_impedance  - Impedance matrix in the modal domain
%     Transformed_admittance - Admittance matrix in the modal domain

	if n_cnds == 3

		%                  Clarke transformation for three phases               %
		% ─────────────────────────────────────────────────────────────────────────
		% Clarke transformation matrix
		TClarke = sqrt(2/3) * [   1         0       1/sqrt(2)
								-1/2    sqrt(3)/2   1/sqrt(2)
								-1/2   -sqrt(3)/2   1/sqrt(2)];

		% Inverse Clarke transformation matrix
		TClarke_inv = sqrt(2/3) * [     1         -1/2        -1/2;
										0       sqrt(3)/2  -sqrt(3)/2;
									1/sqrt(2)   1/sqrt(2)   1/sqrt(2)];

		% Transformation of the impedance matrix
		Transformed_impedance = zeros(n_cnds,n_cnds,N_samples);
		for k = 1:N_samples
			Transformed_impedance(:,:,k) = TClarke_inv * Z_transposed(:,:,k) * TClarke;
		end


		% Transformation of the admittance matrix
		Transformed_admittance = zeros(3,3,N_samples);
		for k = 1:N_samples
			Transformed_admittance(:,:,k) = TClarke_inv * Y_transposed(:,:,k) * TClarke;
		end
		% ─────────────────────────────────────────────────────────────────────────

	elseif n_cnds == 6
		
		%                Clarke transformation for a double circuit             %
		% ─────────────────────────────────────────────────────────────────────────
		% Clarke transformation matrix
		TClarke = [-1/sqrt(6)  1/sqrt(2)   1/sqrt(6)   1/sqrt(6)        0            0
					2/sqrt(6)      0       1/sqrt(6)   1/sqrt(6)        0            0
				   -1/sqrt(6) -1/sqrt(2)   1/sqrt(6)   1/sqrt(6)        0            0
						0          0      -1/sqrt(6)   1/sqrt(6)    1/sqrt(2)    -1/sqrt(6)
						0          0      -1/sqrt(6)   1/sqrt(6)        0         2/sqrt(6)
						0          0      -1/sqrt(6)   1/sqrt(6)   -1/sqrt(2)    -1/sqrt(6)];

		% Inverse Clarke transformation matrix
		TClarke_inv = [-1/sqrt(6)  2/sqrt(6)  -1/sqrt(6)       0            0            0 
						1/sqrt(2)      0      -1/sqrt(2)       0            0            0 
						1/sqrt(6)  1/sqrt(6)   1/sqrt(6)  -1/sqrt(6)   -1/sqrt(6)   -1/sqrt(6);
						1/sqrt(6)  1/sqrt(6)   1/sqrt(6)   1/sqrt(6)    1/sqrt(6)    1/sqrt(6);
						   0           0           0       1/sqrt(2)        0       -1/sqrt(2);
						   0           0           0      -1/sqrt(6)    2/sqrt(6)   -1/sqrt(6)];

		% Transformation of the impedance matrix
		Transformed_impedance = zeros(n_cnds,n_cnds,N_samples);
		for k = 1:N_samples
			Transformed_impedance(:,:,k) = TClarke_inv * Z_transposed(:,:,k) * TClarke;
		end


		% Transformation of the admittance matrix
		Transformed_admittance = zeros(n_cnds,n_cnds,N_samples);
		for k = 1:N_samples
			Transformed_admittance(:,:,k) = TClarke_inv * Y_transposed(:,:,k) * TClarke;
		end
		% ─────────────────────────────────────────────────────────────────────────
	end

end