%Butterworth filter
Wc = 1.0465;        % cut-off frequency
N = 11;             % order of the butterworth filter

%Poles calculation
poles = zeros(1,N);
x = zeros(1,N);
y = zeros(1,N);
for j=0:(N-1)
    x(j+1) = Wc*cos((pi+2*j*pi)/(2*N)+(pi/2));
    y(j+1) = Wc*sin((pi+2*j*pi)/(2*N)+(pi/2));
    poles(j+1) = x(j+1) + 1i*y(j+1);
end

% stopband and passband edges in discrete domain
fp1 = 98;
fs1 = 103;
fs2 = 143;
fp2 = 148;

f_samp = 425;                   % sampling frequency      
wp1 = tan(pi*fp1/f_samp);       % corresponding Analog frequencies
ws1 = tan(pi*fs1/f_samp); 
ws2 = tan(pi*fs2/f_samp);
wp2 = tan(pi*fp2/f_samp);

% parameters of the frequency transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

% evaluating the numerator and denominator polynomials of the Analog Transfer function
[numerator,denominator] = zp2tf([],poles,Wc^N);

syms s z;
Analog_LPF_TF(s) = poly2sym(numerator,s)/poly2sym(denominator,s);% Analog lowpass transfer function
sL = (B*s)/(s*s + W0*W0);
Analog_BSF_TF(s) = Analog_LPF_TF(sL);% Analog bandstop transfer function
Discrete_BSF_TF(z) = Analog_BSF_TF((1-(z^(-1)))/(z^(-1)+1));% discrete bandstop transfer function

[numerator_s, denominator_s] = numden(Analog_BSF_TF(s));
numerator_s = sym2poly(numerator_s);                     % vector of coefficients     
denominator_s = sym2poly(denominator_s);
Den_lead_coeff_s = denominator_s(1);  %Normalizing factor
denominator_s = denominator_s/Den_lead_coeff_s; %normalized
numerator_s = numerator_s/Den_lead_coeff_s;


[numerator_z, denominator_z] = numden(Discrete_BSF_TF(z));                    
numerator_z = sym2poly(numerator_z);        % vector of coefficients
denominator_z = sym2poly(denominator_z);
Den_lead_coeff_z = denominator_z(1);    %Normalizing factor
denominator_z = denominator_z/Den_lead_coeff_z; %normalized
numerator_z = numerator_z/Den_lead_coeff_z;
fvtool(numerator_z,denominator_z)       % magnitude response and phase response with frequency

[H,f] = freqz(numerator_z,denominator_z, 425e3);    % final function
plot(f,abs(H))                  % plotting magnitude response with frequency
grid
