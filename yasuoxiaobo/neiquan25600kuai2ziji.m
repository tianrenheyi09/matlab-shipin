clear all;
close all;
vib=textread('E:\25600数据\REC3953_ch3.txt','','headerlines',16);
pul=textread('E:\25600数据\REC3953_ch2.txt','','headerlines',16);
pluse1=vib(:,1);
vdata=vib(:,2);
pul=pul(:,2);pul=-pul;threshold=abs(1)/3;
fs=25600;
figure;
plot(pluse1,vdata)
xlabel('t/s');ylabel('amplitude/g');title('原始信号');
%%%取降速信号处理：
vdata=vdata(fs*0.5:fs*5.5);%降速7-12
pul=pul(fs*0.5:fs*5.5);
n=length(vdata);
t=(0:n-1)/fs;
% f=(0:n-1)*fs/n;
figure;
plot(t,vdata)
xlabel('t/s');ylabel('amplitude/g');title('截取信号');

sig=vdata;
%%延拓%%
% sig=sig0';
% [sig,startnode,xlength,tnode]=duichengyantuo(sig,fs);
% t=(0:length(sig)-1)/fs;
%% 峭度滤波
% nlevel= 5;
% opt1=1;opt=2;
% figure()
% sig= Fast_Kurtogram2(sig,nlevel,fs,opt1,opt);
%% 低通滤波 
% fp=2000;fst=2100;%滤波器的参数设置
% sig=sig-mean(sig);
% sig=Noiseditong(sig,fs,fp,fst); 

fp=[3100,3700];fst=[2900,3900];%滤波器的参数设置
sig=sig-mean(sig);
sig=Noisedaitong(sig,fs,fp,fst);

% vdata=sig;

figure;
plot(t,sig);
xlabel('t/s');ylabel('amplitude/g')
title('滤波后信号')
sig=sig-mean(sig); 

% vdata=sig;
% siglv=sig';
% %%
% sig=siglv;
fs2=1024;
chongcaiyanlv=floor(fs/fs2);
fs_origin=fs;
pp=1:(floor((length(sig)-1)/chongcaiyanlv)+1);%%新的采样频率下对应的点数位置
 pp=(pp-1)*chongcaiyanlv+1;
sig1=sig(pp);
t2=t(pp);%%%t2
figure;
subplot(211);plot(t,sig);xlabel('t/s');ylabel('amplitude/g');title('滤波后信号');
subplot(212);plot(t2,sig1);xlabel('t/s');ylabel('amplitude/g');title('降采频信号');%%%%t2
set(gcf,'unit','centimeters','position',[3 5 13.5 9])
% xlabel('t/s');ylabel('amplitude');
fs=fs2; 
clear 　chongcaiyanlv  pp fst;clear chongcaiyanlv fs2;pack;
%%%降低采频
%%%出口：变采后的时间fs,t2,sig endfreq fp
%%%出口：原来的fs_origin t vdata 及chishu参数
%************%１）低通滤波及滤波后降频采样**************************************

%%
jieduan1=1.5;jieduan2=3.5;%1-2可以
%*************%２）计算WVD后要截取的信号指针****************************************
%%% 入口：变采后的时间fs,t2,sig  fp
%%% 入口：原来的fs_origin t vdata 及chishu参数
%%% 目的：在信号两端留下１s时长以便清除WVD端点效应
sig=sig1';
%%%%WVD去端点效应再截断两边各1s
% a=abs(t2-jieduan1);
[~,ind1]=min(abs(t2-jieduan1));[~,ind2]=min(abs(t2-jieduan2));
sig=sig-mean(sig);
 sig=sig';
 specnum=length(sig);
 sig=hilbert(sig);
 sig=abs(sig);                        
 sig=sig-mean(sig);
 sig11=hilbert(sig);
g=tftb_window(round(length(sig11)/20+1),'hamming'); h=tftb_window(round(length(sig11)/20+1),'hamming'); 
[tfr,~,f] = tfrspwv(sig11,1:length(sig11),specnum,g,h);

% [tfr,tt,f]=tfrstft(sig11,1:length(sig11),specnum,g,h);%短时傅里叶变换

f=f*fs;

tfr=tfr(:,ind1:ind2);

t1=t2(ind1:ind2);

sig2=sig(ind1:ind2);
% f=f(1:end/2)*fs;%%%%
% tfr=tfr(1:end/2,ind1:ind2);%%tfr(:,ind1:ind2);tfr(1:end/2,ind1:ind2);
% t1=t2(ind1:ind2);

% figure;
% %mesh(t1,f,tfr);
% pcolor(t1,f,tfr);
% ylim([0,250]);
% shading interp;
% colorbar;
% xlabel('t/s');ylabel('f/Hz');
% title('振动信号的SPWV时频分布图')
% set(gcf,'unit','centimeters','position',[3 5 13.5 9])


clear ind1 ind2 fs specnum  tfr_stft; 

%%% 出口：原来的fs_origin t vdata 及chishu参数
%%% 出口：变采后的时间fs,t2,sig  fp
%%% 出口：只是多加了两个截断后时间对齐标志ind1　ind2
%%
%*************%４）vitebi搜索频率路径****************************************
%%%入口参数：时间轴向量t1 频率轴向量f 时频分布tfr  起始频率startfreq 
%%%入口参数：原始振动信号的vdata,fs_origin ,t及齿数参数chishu 
% %%%%%%%划分频率轴格
df=(f(2)-f(1));meiMge=floor(1/df); 
[tfr2,~,ind3]=pinglvfengduan(abs(tfr),f,meiMge);
clear tfr;tfr=tfr2;
clear tfr2 
%%%%%%%划分频率轴格
c=5;sigma=3;efxia=2;
[route,val_route]=viterbi3(tfr,c,sigma);
for p=1:size(route,1),   
    temp=route(p,:);
    temp2=zeros(1,length(temp));
    for time_node=1:length(temp),
        temp2(time_node)=(temp(time_node)-1)*meiMge+ind3(temp(time_node),time_node)-1;
    end
    route(p,:)=temp2;
end
clear temp temp2 time_node p;
 for p=1:size(route,1),   
route(p,:)=f(route(p,:));
end
% f1=route(12,:);
% f1=route(16,:);3942 3
f1=route(5,:);%3

% f1=f1/2;
save('f11.mat','f1');
%%%%对计算出来的估计转速插值回原来的fs
fs=fs_origin;
% f1=interp1(t1,f1,t1(1):1/fs:t1(end));%%%%%
% % f2=interp1(t1,f2,t1(1):1/fs:t1(end));
% f1=f1/2;
% t1=t1(1):1/fs:t1(end);%%%%%%
% [ttt,frequency2]=shijispeed(pul,threshold,fs);
% 
% a=polyfit(t1,f1,3);
% f2=polyval(a,t1,3);
% 
% 
% figure;
% hold on;
% plot(t1,f1,'k');
% plot(t1,f2,'r');
% plot(ttt,frequency2);
% xlabel('t/s');ylabel('amplitude');title('估计转速');
% hold off;
%% 积分求相位函数
fs=1024;
f3=f1-f1(1);
phs = cumtrapz(f3)/fs*2*pi;
phs1=cumtrapz(f1)/fs*2*pi;
tt=(0:length(phs)-1)/fs;
phs = phs(:);
figure;
plot(phs);
%% 广义解调
y = exp(j*(-phs)).*sig2;
n1=length(y);

%% 时频图
% myWindowFT(y,1,500,200,'Hanning',fs,'STFT',1)
[tfr,~,f]=tfrspwv(y,1:n1,n1,hamming(251),hamming(251));
f=f*fs;
% figure;
% contour(tt,f,abs(tfr));
% % ylim([0,250]);
% shading interp;
% colorbar;
% xlabel('t/s');ylabel('f/Hz');
% set(gcf,'unit','centimeters','position',[3 5 13.5 9])
%% 带通滤波
[yy] = mybpf(y,fs,45,54);
n2=length(yy);
[tfr,~,f]=tfrspwv(yy,1:n2,n2,hamming(251),hamming(251));

f=f*fs;
% figure;
% contour(tt,f,abs(tfr));
% shading interp;
% colorbar;
% xlabel('t/s');ylabel('f/Hz');
% set(gcf,'unit','centimeters','position',[3 5 13.5 9])
%% GFT反变换
y3=yy.*exp(-j*(-phs));
[tfr,~,f]=tfrspwv(y3,1:n1,n1,hamming(251),hamming(251));
f=f*fs;
% figure;
% pcolor(tt,f,abs(tfr));
% shading interp;
% colorbar;
% xlabel('t/s');ylabel('f/Hz');
% set(gcf,'unit','centimeters','position',[3 5 13.5 9])


%% 取GFT反变换后的实部
y4=real(y3);
figure;
plot(tt,y4);
xlabel('t/s');ylabel('amplitude/g');
set(gcf,'unit','centimeters','position',[3 5 13.5 9])
%% 求瞬时相位
y5=hilbert(y4);
yh1=unwrap(angle(y5));%2次谐波瞬时相位
% yh1=yh1/2;%%轴瞬时相位
figure;
plot(tt,yh1,'r','linewidth',2);
% hold on
% plot(tt,phs1);
xlabel('t/s');ylabel('phase/rad');
set(gcf,'unit','centimeters','position',[3 5 13.5 9])

yh1=yh1/2;%%轴瞬时相位
phs1=phs1/2;
f1=f1/2;
yhd=fs*diff(yh1)/(2*pi);%瞬时频率
figure;
plot(yhd);



%% 重采样
[~,ind1]=min(abs(t-jieduan1));[~,ind2]=min(abs(t-jieduan2)); 
vdata1=vdata(ind1:ind2);
x1=vdata1;
% x1=hilbert(x1);

fs=25600;
% f1=interp1(t1,f1,t1(1):1/fs:t1(end));%%%%%
% f2=interp1(t1,f2,t1(1):1/fs:t1(end));
phs1=interp1(tt,phs1,tt(1):1/fs:tt(end));
yh1=interp1(tt,yh1,tt(1):1/fs:tt(end));

% fs=25600;
data=2*pi*min(f1)/fs;
tangle=0:data:phs1(end);
tangle1=0:data:yh1(end);

array_angle_amp=interp1(phs1,x1,tangle,'spline');%插值
array_angle_amp1=interp1(phs1,x1,tangle1,'spline');%插值
figure;
plot(tangle,array_angle_amp);%重采样信号
% set(gcf,'unit','centimeters','position',[3 5 13.5 9])
% figure;
% plot(tangle1,array_angle_amp1);
% xlabel('angle/rad');ylabel('amplitude/g');

array_angle_amp1=array_angle_amp1(1:end-0.1*fs);%消除波动
% array_angle_amp1=hilbert(array_angle_amp1);
% array_angle_amp1=abs(array_angle_amp1);
% array_angle_amp1=array_angle_amp1-mean(array_angle_amp1);
figure;
plot(tangle1(1:end-0.1*fs),array_angle_amp1);
% xlabel('angle/rad');ylabel('amplitude/g');
% set(gcf,'unit','centimeters','position',[3 5 13.5 9])


array_angle_amp=hilbert(array_angle_amp1);
array_angle_amp=abs(array_angle_amp);
array_angle_amp=array_angle_amp-mean(array_angle_amp);

angle_dom_ffty=abs(fft(array_angle_amp))*2/length(array_angle_amp);
 delt_order=2*pi/(length(angle_dom_ffty)*data);
angle_dom_fx=(0:length(angle_dom_ffty)-1)*delt_order;%类似时域的（0：n-1)*n/fs
freq=angle_dom_fx(1:length(angle_dom_fx)/2);%阶次
ayy=angle_dom_ffty(1:length(angle_dom_fx)/2);
figure;
plot(freq,ayy);
xlabel('order');ylabel('amplitude/g');
xlim([0,20]);
set(gcf,'unit','centimeters','position',[3 5 13.5 9])

ffs=2*pi/data;
[freq,ayy,vec_num]=zijijieci(array_angle_amp1,ffs,data);
figure;
plot(freq,ayy);
xlabel('order');ylabel('amplitude/g');
xlim([0,20]);
set(gcf,'unit','centimeters','position',[3 5 13.5 9])

