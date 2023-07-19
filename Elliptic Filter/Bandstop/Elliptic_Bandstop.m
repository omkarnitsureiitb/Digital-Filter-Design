% Elliptic Bandstop Filter
Delta = 0.15;

% Specifications Borrowed From The Butterworth Bandstop Filter Design
W0 = 1.30889;
B = 1.0511;
Wp = 1;
Ws = 1.2447;

% Now we find the minimum allowed passband Magnitude and maximum stopband
% magnitude
A_Pb_min = 1 - Delta;
A_Sb_max = Delta;

% Calculating the ripple parameters epsilon_Pb, epsilon_Sb
epsilon_Pb = sqrt((1/(A_Pb_min)^2) - 1);
epsilon_Sb = sqrt(1/((A_Sb_max)^2) - 1);

% Calculating the Selectivity(k) and Discrimination(k1) parameters
k = Wp/Ws;
k1 = epsilon_Pb/epsilon_Sb;

% Solving for the Elliptic Integrals K, K', K1, K1' using the Matlab
% Function ellipk
[K, Kp] = ellipk(k);
[K1, K1p] = ellipk(k1);

% Solving for the Order of Elliptic Filter
N_float = (K/K1)*(K1p/Kp);

% We get N_float = 2.7341 and thus taking the minimum integer greater than N, we
% get Order of Elliptic Filter, N = 3
N = 3;

% But now N = 3 does not satisfy the Degree equation, so we recalculate the
% Selectivity Parameter(k) by putting new N and the original k1
k = ellipdeg(N, k1);

% We can write N = 2L + r, where L is the number of Quadratic Factors and r
% is the number of real poles
% We have N = 4 and thus L = 1, r = 1
L = 1;
r = 1;

% We now calculate the matrix of Scaling Factor u of z in the Jacobian
% Elliptic Integrals
u = (2*1 - 1)/N;

% We now calculate zeros of Fn(w) in the Elliptic Filter Magnitude squared
% equation
zeta_i = cde(u,k);

% Now we evaluate the zeros of the Elliptic Filter in the Stopband unlike
% the Butterworth or Chebyschev Filter where there were no zeros
Zeroes = 1i*(Wp/(k*zeta_i));

% Calculating the parameter v0 to evaluate the poles of the Elliptic Filter
v0 = -1i*asne(1i/epsilon_Pb, k1)/N;

% Now evaluating the poles of the Transfer Function
% Also as N is odd we have a real pole
Poles = Wp*1i*cde(u - 1i*v0, k);
Pole_real = Wp * 1i*sne(1i*v0, k);

% As there are 4 zeros and 4 poles we have 4th Degree polynomials in both
% Numerator and Denominator of the Transfer Function
num_coeff = zeros(L, 3);
den_coeff = zeros(L, 3);
num_coeff(:, 3) = 1;
den_coeff(:, 3) = 1;

num_coeff(2) = -2*real(1/Zeroes);
num_coeff(1) = abs(1/Zeroes)^2;
den_coeff(2) = -2*real(1/Poles);
den_coeff(1) = abs(1/Poles)^2;

% Here as r = 1, we will have to take care of a Linear term in the
% Denominator of the Transfer Function
den_coeff = [0, -real(1/Pole_real), 1; den_coeff];

% Generating Coefficients of the Denominator and Numerator of the
% Analog_Lowpass Filter Specifications
den1 = den_coeff(1,:);
den2 = den_coeff(2,:);

den_coeff = conv(den1,den2);

% Writing down the Analog Lowpass Transfer Function and then using
% Frequency transformation Converting it into the Bandstop Filter Transfer
% Function
syms s z;
Analog_LPF_TF(s) = poly2sym(num_coeff, s)/poly2sym(den_coeff, s);
sL = (B*s)/(s*s + W0*W0);
Analog_BSF_TF(s) = Analog_LPF_TF(sL);

% Converting the Analog Bandstop Filter to Digital Bandstop Filter using
% Frequency Transformation
s_trans = (1 - z^(-1))/(1 + z^(-1));
Discrete_BSF_TF(z) = Analog_BSF_TF(s_trans);

% Obtaining coefficients of Analog BSF
[num_s, den_s] = numden(Analog_BSF_TF(s));      % Extract Numerator and Denominator of the Transfer Function
num_s = sym2poly(num_s);        % Extracts the coefficients of a polynomial in matrix form
den_s = sym2poly(den_s);

% Dividing Numerator and Denominator by the coefficient of Highest Degree
% of Denominator
norm_coeff = den_s(1);
num_s = num_s/norm_coeff;
den_s = den_s/norm_coeff;

% Obtaining coefficients of Final Digital BSF
[num_z, den_z] = numden(Discrete_BSF_TF(z));    % Extract Numerator and Denominator of the Transfer Function
num_z = sym2poly(num_z);        % Extracts the coefficients of a polynomial in matrix form
den_z = sym2poly(den_z);

% Dividing Numerator and Denominator by the coefficient of Highest Degree
% of Denominator
norm_coeff = den_z(1);
num_z = num_z/norm_coeff;
den_z = den_z/norm_coeff;

% Displays the Final Magnitude and Phase Response of the Digital Bandstop
% Filter designed using the Elliptic Approximation
fvtool(num_z, den_z)

% Magnitude plot in terms of Frequency in kHz
[Hdigital, f] = freqz(num_z, den_z, 1024*1024, 425e3);
plot(f, abs(Hdigital))
grid