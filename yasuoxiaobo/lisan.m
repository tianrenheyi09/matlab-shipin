clc;clear all;close all;

fs=1000;
t=0:1/fs:1;
% sig=cos(2*pi*(10*t.^2+20*t))+cos(2*pi*(50*t.^2+100*t)); %%非平稳
sig=cos(2*pi*(20*t))+cos(2*pi*(100*t)); 
figure;
plot(t,sig);


%% % % 连续小波
%%%%%%%%%%%%%
figure;
cw1=cwt(sig,1:32,'db4','plot');
xlabel('原始信号');
[cw1,sc]=cwt(sig,1:32,'db4','scal');
title('连续变换，绝对系数');
xlabel('时间');ylabel('尺度');
%%%%%%%%%%%

s=sig;
wavename='cmor3-3';
% wavename='morlet';
totalscal=128;
Fc=centfrq(wavename); % 小波的中心频率
c=2*Fc*totalscal;
scals=c./(1:totalscal);
f=scal2frq(scals,wavename,1/fs); % 将尺度转换为频率
coefs=cwt(s,scals,wavename); % 求连续小波系数
figure;
f=(0:length(sig)/2-1)*fs/length(sig);
imagesc(t,f,abs(coefs));
set(gca,'YDir','normal')
colorbar;
xlabel('t/s');
ylabel('f/Hz');
title('小波时频图');


%% %% 离散小波变换
[ca,cd]=dwt(sig,'haar');
figure;
subplot(211);plot(ca);
subplot(212);plot(cd);
%%%%%%%%%%%%%%
figure;
[c,l]=wavedec(s,5,'db10');
%%低频部分
for i=1:5
    qq=wrcoef('a',c,l,'db10',6-i);
    subplot(5,1,i);
    plot(qq);
    ylabel(['a',num2str(6-i)]);
end
%%高频部分
figure;
for i=1:5
    qq=wrcoef('d',c,l,'db10',6-i);
    subplot(5,1,i);
    plot(qq);
    ylabel(['a',num2str(6-i)]);
end
    

    