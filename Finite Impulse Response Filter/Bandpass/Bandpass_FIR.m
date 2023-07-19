% FIR Bandpass Filter Parameters
fsampling = 600000;
Delta = 0.15;

%Bandpass Filter frequency speifications at the two edges of the Passband
fs1 = 106000;
fp1 = 111000;
fp2 = 186000;
fs2 = 191000;

% cutoff Frequencies Wc1, Wc2 for ideal Bandpass Frequency Response
fc1 = (fs1 + fp1)/2;
fc2 = (fs2 + fp2)/2;

Wc1 = (fc1/fsampling)*2*pi;
Wc2 = (fc2/fsampling)*2*pi;

% wT = Ws - Wp
wT = 2*pi*(fp1 - fs1)/fsampling;

% Calculating N according to the following formula
% N >= (A - 8)/(2*2.285*wT)
A = -20*log10(Delta);
N_float = (A - 8)/(2*2.285*wT);

% We get N_float = 35.4313, thus we take N = 36
N = 36;

% Calculating Alpha for the given tolerance level
if(A < 21)
    Alpha = 0;
elseif((A >= 21) && (A <= 50))
    Alpha = 0.5842*(A - 21)^0.4 + 0.07886*(A - 21);
else
    Alpha = 0.1102*(A - 8.7);
end

Beta = Alpha/N;

% Now let us find the ideal Bandpass impulse response by using the fact 
% that the lowpass filter is an LTI system so we can subtract 1 Lowpass
% frequency response from another to get Bandpass
Lowpass_ideal_Wc2 = zeros(1,2*N + 1);
Lowpass_ideal_Wc1 = zeros(1,2*N + 1);
for index = -N:N
    if(index ~= 0)
        Lowpass_ideal_Wc2(index + (N + 1)) = sin(Wc2*index)/(pi*index);
        Lowpass_ideal_Wc1(index + (N + 1)) = sin(Wc1*index)/(pi*index);
    else
        Lowpass_ideal_Wc2(N + 1) = Wc2/pi;
        Lowpass_ideal_Wc1(N + 1) = Wc1/pi;
    end
end
Bandpass_ideal = Lowpass_ideal_Wc2 - Lowpass_ideal_Wc1;

%Kaiser Window of length 2N + 1 with shape parameter beta calculated above
kaiser_win_coeff = kaiser((2*N + 1),Beta);

% Calculation the Bandpass FIR filter coefficients
Bandpass_FIR_Filter = zeros(1, (2*N + 1));
for index = 1:(2*N + 1)
    Bandpass_FIR_Filter(index) = Bandpass_ideal(index)*kaiser_win_coeff(index);
end

%Evaluating frequency response of the Filter
fvtool(Bandpass_FIR_Filter);

%magnitude response in deciBels
[H,f] = freqz(Bandpass_FIR_Filter,1,1024, fsampling);
plot(f,abs(H))
grid