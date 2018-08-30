clc
clear all
close all
% [x, Fs] = wavread('Hum.wav');
% Ts = 1/Fs;
% x = x(1:6000);
Ts = 0.001;
Fs = 1/Ts;
t=0:Ts:1;
x = sin(2*pi*10*t) + sin(2*pi*50*t) + sin(2*pi*100*t) + 0.1*randn(1, length(t));
imf = emd(x);
plot_hht(x,imf,1/Fs);
k = 4;
y = imf{k};
N = length(y);
t = 0:Ts:Ts*(N-1);
[yenvelope, yfreq, yh, yangle] = HilbertAnalysis(y, 1/Fs);
yModulate = y./yenvelope;
[YMf, f] = FFTAnalysis(yModulate, Ts);
Yf = FFTAnalysis(y, Ts);
figure
subplot(321)
plot(t, y)
title(sprintf('IMF%d', k))
xlabel('Time/s')
ylabel(sprintf('IMF%d', k));
subplot(322)
plot(f, Yf)
title(sprintf('IMF%d的频谱', k))
xlabel('f/Hz')
ylabel('|IMF(f)|');
subplot(323)
plot(t, yenvelope)
title(sprintf('IMF%d的包络', k))
xlabel('Time/s')
ylabel('envelope');
subplot(324)
plot(t(1:end-1), yfreq)
title(sprintf('IMF%d的瞬时频率', k))
xlabel('Time/s')
ylabel('Frequency/Hz');
subplot(325)
plot(t, yModulate)
title(sprintf('IMF%d的调制信号', k))
xlabel('Time/s')
ylabel('modulation');
subplot(326)
plot(f, YMf)
title(sprintf('IMF%d调制信号的频谱', k))
xlabel('f/Hz')
ylabel('|YMf(f)|');

figure;
plot_hht(x,imf,Ts)