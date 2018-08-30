% clc;clear;close all
% t1=0:1/40:7.5-1/40;
% t2=7.5:1/40:25.6-1/40;
% t=[t1,t2];
% x1=sin(1.6*2*pi*t1)+sin(3*2*pi*t1+0.04*2*pi*t1.^2+4*sin(0.1*2*pi*t1));
% x2=sin(1.6*2*pi*(t2-7.5)+0.055*2*pi*(t2-7.5).^2)+sin(3*2*pi*t2+0.04*2*pi*t2.^2+4*sin(0.1*2*pi*t2));
% x=[x1,x2];
% N=length(x);
% y=awgn(x,11);
% fs=40;fc=3;totalscal=200;
% FreqBins=linspace(1/fs,6/fs,totalscal);
% fw=FreqBins*fs;
% scales=fc./FreqBins;


clear;clc;close all
SamFreq=2048;
N=3*SamFreq;
fs=SamFreq;
t=0:1/fs:(N-1)/fs;

s_f1=(2.5*t.^2+10*t);s_f2=(12.5*t.^2+50*t);
f1=5*t+20;f2=50+30*t;
f1_zheng=f1;f2_zheng=f2;
figure(1);
plot(t,f1,'r');hold on;plot(t,f2,'b');
sig=cos(2*pi*s_f1)+sin(2*pi*s_f2); %另别的sig


%%%降低采频
fs2=1024;chongcaiyanlv=floor(fs/fs2);
fs_origin=fs;
pp=1:(floor((length(sig)-1)/chongcaiyanlv)+1);%降低采频后信号的长度
 pp=(pp-1)*chongcaiyanlv+1;%每隔chonggcaiyanlv间隔取一个点，长度还是length(pp)
sig=sig(pp);
t2=t(pp);
figure;plot(t2,sig);title('滤波后再降频采样信号');xlabel('t/s');ylabel('幅值');
fs=fs2; 
clear fs2　chongcaiyanlv  pp fst;clear chongcaiyanlv fs2;pack;
%%%降低采频
%%%出口：变采后的时间fs,t2,sig endfreq fp
%%%出口：原来的fs_origin t vdata 及chishu参数
%************%１）低通滤波及滤波后降频采样**************************************


s=sig;
wavename='cmor3-3';
totalscal=512;
Fc=centfrq(wavename); % 小波的中心频率
c=2*Fc*totalscal;
scals=c./(1:totalscal);
f=scal2frq(scals,wavename,1/fs); % 将尺度转换为频率

figure;
scales=scals;
coes_y=cwt(s,scales,wavename);
imagesc(t2,scales,abs(coes_y));
figure;
imagesc(t2,f,abs(coes_y));
set(gca,'YDir','normal')
colorbar;



%% 局部极大值搜索

% f=(0:length(sig)/2-1)*fs/length(sig);

fw=f;
N=length(s);
a=zeros(1,N);
b=zeros(1,N);
b_p=zeros(1,N);
[a(1),b(1)]=max(abs(coes_y(1:end,1)));
temp1=20;%%开始搜索的频率
b_p(1)=b(1)+temp1;
for i=2:N
    [a(i),b(i)]=max(abs(coes_y(temp1+b(i-1)-35:temp1+b(i-1)+35,i)));
    temp1=b_p(i-1)-35-1;
    b_p(i)=b(i)+temp1;
   coes_y(temp1+b(i-1)-10:temp1+b(i-1)+10,i)=0;
end

ff1=fw(b_p);
figure;
hold on;
% plot(t2,f2,'r');
% p=polyfit(t2,ff2,1);
% y0=polyval(p,t2);
plot(t,f2,'r',t2,ff1);




[a(1),b(1)]=max(abs(coes_y(1:end,1)));
temp1=10;%%开始搜索的频率
b_p(1)=b(1)+temp1;
for i=2:N
    [a(i),b(i)]=max(abs(coes_y(temp1+b(i-1)-5:temp1+b(i-1)+5,i)));
    temp1=b_p(i-1)-5-1;
    b_p(i)=b(i)+temp1;
   
end

ff2=fw(b_p);
figure;
hold on;
% plot(t2,f2,'r');
% p=polyfit(t2,ff2,1);
% y0=polyval(p,t2);
plot(t,f2,'r',t2,ff2);

figure;
plot(t,f1,'r',t,f2,'r');
hold on;
plot(t2,ff1,t2,ff2);