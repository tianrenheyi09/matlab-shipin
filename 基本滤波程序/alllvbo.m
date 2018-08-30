% function y= myditong(x,f1,f3,fs,rp,rs)
clc;clear all;close all;
fs=2000;
t=0:1/fs:0.3;
y=sin(2*pi*300*t)+sin(2*pi*100*t);
figure;
plot(t,y);
%% 巴特沃斯——低通滤波器设计
% rp=0.1;rs=50;
% Wp = 2*120/fs; Ws = 2*200/fs;
% [n,Wn] = buttord(Wp,Ws,rp,rs);
% [b,a] = butter(n,Wn);
% figure;
% freqz(b,a,512,fs); 
% y1=filter(b,a,y);
% figure;
% plot(t,y1);
% 
% n=length(y1);
% f=(0:n/2-1)*fs/n;
% s=abs(fft(y1-mean(y1)));
% figure;
% plot(f,s(1:n/2));

%% 巴特沃斯——高通滤波器设计
% rp=0.1;rs=50;
% Wp = 2*200/fs; Ws = 2*250/fs;
% [n,Wn] = buttord(Wp,Ws,rp,rs);
% [b,a] = butter(n,Wn,'high');
% figure;
% freqz(b,a,512,fs); 
% y1=filter(b,a,y);
% figure;
% plot(t,y1);
% 
% n=length(y1);
% f=(0:n/2-1)*fs/n;
% s=abs(fft(y1-mean(y1)));
% figure;
% plot(f,s(1:n/2));

%% 巴特沃斯——带通滤波器设计

% Wp = [250 350]/(fs/2); Ws = [200 380]/(fs/2);
% rp = 0.1; rs = 40;
% [n,Wn] = buttord(Wp,Ws,rp,rs);
% % Returns n =16; Wn =[0.1198 0.4005];
% [b,a] = butter(n,Wn);
% freqz(b,a,512,fs)
% 
% y1=filter(b,a,y);
% figure;
% plot(t,y1);
% 
% n=length(y1);
% f=(0:n/2-1)*fs/n;
% s=abs(fft(y1-mean(y1)));
% figure;
% plot(f,s(1:n/2));


%% 切比雪夫低通
% Wp = 2*120/fs; Ws = 2*200/fs;
% rp = 3; rs = 60;
% [n,Wp] = cheb1ord(Wp,Ws,rp,rs)
% % Returns n = 4 Wp =0.0800
% [b,a] = cheby1(n,rp,Wp);
% figure;
% freqz(b,a,512,fs);
% y1=filter(b,a,y);
% figure;
% plot(t,y1);
% n=length(y1);
% f=(0:n/2-1)*fs/n;
% s=abs(fft(y1-mean(y1)));
% figure;
% plot(f,s(1:n/2));

%% 切比雪夫带通

% Wp = [250 350]/(fs/2); Ws = [200 380]/(fs/2);
% rp = 3; rs = 40;
% [n,Wp] = cheb1ord(Wp,Ws,rp,rs)
% % Returns n =7 Wp =[0.1200    0.4000]
% [b,a] = cheby1(n,rp,Wp);
% freqz(b,a,512,fs);
% 
%  y1=filter(b,a,y);
% figure;
% plot(t,y1);
% n=length(y1);
% f=(0:n/2-1)*fs/n;
% s=abs(fft(y1-mean(y1)));
% figure;
% plot(f,s(1:n/2));

%% 切比雪夫高通
% Wp = 2*200/fs; Ws = 2*250/fs;
% rp = 3; rs = 60;
% [n,Wp] = cheb1ord(Wp,Ws,rp,rs)
% % Returns n = 4 Wp =0.0800
% [b,a] = cheby1(n,rp,Wp,'high');
% figure;
% freqz(b,a,512,fs);
% y1=filter(b,a,y);
% figure;
% plot(t,y1);
% n=length(y1);
% f=(0:n/2-1)*fs/n;
% s=abs(fft(y1-mean(y1)));
% figure;
% plot(f,s(1:n/2));

%% FIR低通
fp=120;fst=150;
wp=2*fp*pi/fs;ws=2*fst*pi/fs;
deltaw=abs(ws-wp);
N=ceil(6.6*pi/deltaw);
N=N+mod(N,2);
wind=(hamming(N+1))';
Wn=(fp+fst)/fs;
% b=fir1(N,Wn,wind,'high');
b=fir1(N,Wn);%% b=fir1(50,Wn,'high');
figure;
freqz(b,1,512,fs);


y1=filter(b,1,y);
figure;
plot(t,y1);

n=length(y1);
f=(0:n/2-1)*fs/n;
s=abs(fft(y1-mean(y1)));
figure;
plot(f,s(1:n/2));
% % 
%% fir高通
fp=200;fst=250;
wp=2*fp*pi/fs;ws=2*fst*pi/fs;
deltaw=abs(ws-wp);
N=ceil(6.6*pi/deltaw);
N=N+mod(N,2);
wind=(hamming(N+1))';
Wn=(fp+fst)/fs;
% b=fir1(N,Wn,wind,'high');
b=fir1(N,Wn,'high');%% b=fir1(50,Wn,'high');
figure;
freqz(b,1,512,fs);


y1=filter(b,1,y);
figure;
plot(t,y1);

n=length(y1);
f=(0:n/2-1)*fs/n;
s=abs(fft(y1-mean(y1)));
figure;
plot(f,s(1:n/2));

%% FIR带通
fp1=250;fp2=320;
fst1=200;fst2=350;
ws1=fst1*pi*2/fs;wp1=2*fp1*pi/fs;
ws2=fst2*pi*2/fs;wp2=2*fp2*pi/fs;
trw=min((wp1-ws1),(ws2-wp2));
M=ceil(6.2*pi/trw);
M=M+mod(M+1,2);
wc1=(ws1+wp1)/2;wc2=(wp2+ws2)/2;
fc1=wc1/pi;fc2=wc2/pi;
h1=fir1(M-1,[fc1,fc2],hanning(M)');

figure;
freqz(h1,1,512,fs);

y1=filter(h1,1,y);
figure;
plot(t,y1);
n=length(y1);
f=(0:n/2-1)*fs/n;
s=abs(fft(y1-mean(y1)));
figure;
plot(f,s(1:n/2));




