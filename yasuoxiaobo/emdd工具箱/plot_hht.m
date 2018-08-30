function plot_hht(x,imf,Ts)
% Plot the HHT.
% :: Syntax
%    The array x is the input signal and Ts is the sampling period.
%    Example on use: [x,Fs] = wavread('Hum.wav');
%                    plot_hht(x(1:6000),1/Fs);
% Func : emd
% imf = emd(x);
for k = 1:length(imf)
    b(k) = sum(imf{k}.*imf{k});
    th   = unwrap(angle(hilbert(imf{k})));  % 相位
    d{k} = diff(th)/Ts/(2*pi);          % 瞬时频率
end
[u,v] = sort(-b);
b     = 1-b/max(b);                     % 后面绘图的亮度控制
N = length(x);
c = linspace(0,(N-2)*Ts,N-1);           % 0:Ts:Ts*(N-2)
for k = v(1:2)                          % 显示能量最大的两个IMF的瞬时频率
    figure
    plot(c,d{k});
    xlim([0 c(end)]);
    ylim([0 1/2/Ts]);
    xlabel('Time/s')
    ylabel('Frequency/Hz');
    title(sprintf('IMF%d', k))
end
% 显示各IMF
M = length(imf);
N = length(x);
c = linspace(0,(N-1)*Ts,N);             % 0:Ts:Ts*(N-1)
for k1 = 0:4:M-1
    figure
    for k2 = 1:min(4,M-k1)
        subplot(4,2,2*k2-1)
        plot(c,imf{k1+k2})
        set(gca,'FontSize',8,'XLim',[0 c(end)]);
        title(sprintf('第%d个IMF', k1+k2))
        xlabel('Time/s')
        ylabel(sprintf('IMF%d', k1+k2));
       
        subplot(4,2,2*k2)
        [yf, f] = FFTAnalysis(imf{k1+k2}, Ts);
        plot(f, yf)
        title(sprintf('第%d个IMF的频谱', k1+k2))
        xlabel('f/Hz')
        ylabel('|IMF(f)|');
    end
 end
figure
subplot(211)
plot(c,x)
set(gca,'FontSize',8,'XLim',[0 c(end)]);
title('原始信号')
xlabel('Time/s')
ylabel('Origin');
subplot(212)
[Yf, f] = FFTAnalysis(x, Ts);
plot(f, Yf)
title('原始信号的频谱')
xlabel('f/Hz')
ylabel('|Y(f)|');