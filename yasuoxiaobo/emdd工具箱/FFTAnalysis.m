% ÆµÆ×·ÖÎö
function [Y, f] = FFTAnalysis(y, Ts)
Fs = 1/Ts;
L = length(y);
NFFT = 2^nextpow2(L);
y = y - mean(y);
Y = fft(y, NFFT)/L;
Y = 2*abs(Y(1:NFFT/2+1));
f = Fs/2*linspace(0, 1, NFFT/2+1);
end