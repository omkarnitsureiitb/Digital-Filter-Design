% Chebyschev Filter

% Necessary Parameters
Delta = 0.15;
D1 = 1/(Delta^2) - 1;
epsilon = sqrt(D1);
N = 6;

% Poles of the Analog Lowpass Filter using WolframAlpha
p1 = -0.05458 + 1i*0.98717;
p2 = -0.05458 - 1i*0.98717;
p3 = -0.2037 + 1i*0.26451;
p4 = -0.2037 - 1i*0.26451;
p5 = -0.14912 + 1i*0.72266;
p6 = -0.14912 - 1i*0.72266;

poles = [p1, p2, p3, p4, p5, p6];

% Sampling Frequency in kHz
f_samp = 600;

% Digital Filter Band Edge Specifications in kHz
fs1 = 106;
fp1 = 111;
fp2 = 186;
fs2 = 191;

% Normalized Frequencies
fs1n = (fs1/f_samp)*(2*pi);
fs2n = (fs2/f_samp)*(2*pi);
fp1n = (fp1/f_samp)*(2*pi);
fp2n = (fp2/f_samp)*(2*pi);

% Using Bilinear Transformation to get the corresponding Analog Specifications
ws1 = tan(fs1n/2);
ws2 = tan(fs2n/2);
wp1 = tan(fp1n/2);
wp2 = tan(fp2n/2);

% DC gain value for the Analog Lowpass Filter Transfer Function
Gain = real(p1*p2*p3*p4*p5*p6)*sqrt(1/(1 + epsilon*epsilon));

% Finding the Coefficients of the numerator and Denominator of the Analog
% Lowpass Transfer function using its poles and zeroes
[num_coeff, den_coeff] = zp2tf([], poles, Gain);

% Parameters for Bandpass to Lowpass Transformation such that wp1, wp2 
% mapped to -+ 1 respectively
W0 = sqrt(wp1*wp2);
B = wp2 - wp1;

% Using the poles of Analog LPF and Frequency Transformation to get the 
% Analog Bandpass Transfer Function
syms s z;
Analog_LPF_TF(s) = poly2sym(num_coeff, s)/poly2sym(den_coeff, s);
sL = (s*s + W0*W0)/(B*s);
Analog_BPF_TF(s) = Analog_LPF_TF(sL);

% Converting the Analog Bandpass Filter to Digital Bandpass Filter using
% Frequency Transformation
s_trans = (1 - z^(-1))/(1 + z^(-1));
Discrete_BPF_TF(z) = Analog_BPF_TF(s_trans);

% Obtaining coefficients of Analog BPF
[num_s, den_s] = numden(Analog_BPF_TF(s));      % Extract Numerator and Denominator of the Transfer Function
num_s = sym2poly(num_s);        % Extracts the coefficients of a polynomial in matrix form
den_s = sym2poly(den_s);

% Dividing Numerator and Denominator by the coefficient of Highest Degree
% of Denominator
norm_coeff = den_s(1);
num_s = num_s/norm_coeff;
den_s = den_s/norm_coeff;

% Obtaining coefficients of Final Digital BPF
[num_z, den_z] = numden(Discrete_BPF_TF(z));    % Extract Numerator and Denominator of the Transfer Function
num_z = sym2poly(num_z);        % Extracts the coefficients of a polynomial in matrix form
den_z = sym2poly(den_z);

% Dividing Numerator and Denominator by the coefficient of Highest Degree
% of Denominator
norm_coeff = den_z(1);
num_z = num_z/norm_coeff;
den_z = den_z/norm_coeff;

% Displays the Final Magnitude and Phase Response of the Digital Bandpass
% Filter designed using the Chebyschev Approximation
fvtool(num_z, den_z)

% Magnitude plot in terms of Frequency in kHz
[Hdigital, f] = freqz(num_z, den_z, 1024*1024, 600e3);
plot(f, abs(Hdigital))
grid