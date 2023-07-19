% FIR Bandstop Filter Parameters
fsampling = 425000;
Delta = 0.15;

%Bandstop Filter frequency speifications at the two edges of the Stopband
fp1 = 98000;
fs1 = 103000;
fs2 = 143000;
fp2 = 148000;

% cutoff Frequencies Wc1, Wc2 for ideal Bandstop Frequency Response
fc1 = (fs1 + fp1)/2;
fc2 = (fs2 + fp2)/2;

Wc1 = (fc1/fsampling)*2*pi;
Wc2 = (fc2/fsampling)*2*pi;

% wT = Ws - Wp
wT = 2*pi*(fs1 - fp1)/fsampling;

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

% Now let us find the ideal Bandstop impulse response by using the fact 
% that the lowpass filter is an LTI system so we can subtract 1 Lowpass
% frequency response from an All-pass system to get the Bandstop impulse
% Response
Lowpass_ideal_Wc2 = zeros(1,2*N + 1);
Lowpass_ideal_Wc1 = zeros(1,2*N + 1);
Lowpass_ideal_pi = zeros(1,2*N + 1);
for index = -N:N
    if(index ~= 0)
        Lowpass_ideal_Wc2(index + (N + 1)) = sin(Wc2*index)/(pi*index);
        Lowpass_ideal_Wc1(index + (N + 1)) = sin(Wc1*index)/(pi*index);
        Lowpass_ideal_pi(index + (N + 1)) = sin(pi*index)/(pi*index);
    else
        Lowpass_ideal_Wc2(N + 1) = Wc2/pi;
        Lowpass_ideal_Wc1(N + 1) = Wc1/pi;
        Lowpass_ideal_pi(N + 1) = 1;
    end
end
Bandstop_ideal = Lowpass_ideal_pi - (Lowpass_ideal_Wc2 - Lowpass_ideal_Wc1);

%Kaiser Window of length 2N + 1 with shape parameter beta calculated above
kaiser_win_coeff = kaiser((2*N + 1),Beta);

% Calculation the Bandpass FIR filter coefficients
Bandstop_FIR_Filter = zeros(1, (2*N + 1));
for index = 1:(2*N + 1)
    Bandstop_FIR_Filter(index) = Bandstop_ideal(index)*kaiser_win_coeff(index);
end

%Evaluating frequency response of the Filter
fvtool(Bandstop_FIR_Filter);

%magnitude response in deciBels
[H,f] = freqz(Bandstop_FIR_Filter,1, 1024, fsampling);
plot(f,abs(H))
grid