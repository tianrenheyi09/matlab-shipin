function [f,z]=fftfenxi(t,y)
L=length(t);N=2^nextpow2(L);
%fft默认计算的信号是从0开始的
t=linspace(t(1),t(L),N);deta=t(2)-t(1);
m=0:N-1;
f=1./(N*deta)*m;
%下面计算的Y就是x(t)的傅里叶变换数值
%Y=exp(i*4*pi*f).*fft(y)%将计算出来的频谱乘以exp(i*4*pi*f)得到频移后[-2,2]之间的频谱值
Y=fft(y);
z=sqrt(Y.*conj(Y));
