clear all;
close all;
vib=textread('E:\code_for1\REC3950_ch3.txt','','headerlines',16);
pluse1=vib(:,1);
vdata=vib(:,2);
fs=25600;
figure;
plot(pluse1,vdata)
xlabel('t/s');ylabel('amplitude');title('原始信号');
%%%取降速信号处理：
vdata=vdata(fs*11:fs*13);
n=length(vdata);
t=(0:n-1)/fs;
% f=(0:n-1)*fs/n;
figure;
plot(t,vdata)
xlabel('t/s');ylabel('amplitude');title('截取信号');

sig=vdata;
%%延拓%%
% sig=sig0';
% [sig,startnode,xlength,tnode]=duichengyantuo(sig,fs);
% t=(0:length(sig)-1)/fs;
%% 峭度滤波
% nlevel= 2;
% opt1=2;opt=1;
% figure()
% sig= Fast_Kurtogram2(sig,nlevel,fs,opt1,opt);
%% 低通滤波 
fp=1000;fst=3000;%滤波器的参数设置
sig=sig-mean(sig);
sig=Noiseditong(sig,fs,fp,fst); 
figure;
plot(t,sig);
xlabel('t/s');ylabel('amplitude')
title('滤波后信号')
sig=sig-mean(sig); 
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
plot(t2,sig1);%%%%t2
title('降采频信号');
xlabel('t/s');ylabel('amplitude');
fs=fs2; 
clear 　chongcaiyanlv  pp fst;clear chongcaiyanlv fs2;pack;
%%%降低采频
%%%出口：变采后的时间fs,t2,sig endfreq fp
%%%出口：原来的fs_origin t vdata 及chishu参数
%************%１）低通滤波及滤波后降频采样**************************************

%%
jieduan1=1;jieduan2=2;
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


f=f*fs;

tfr=tfr(:,ind1:ind2);

t1=t2(ind1:ind2);


% f=f(1:end/2)*fs;%%%%
% tfr=tfr(1:end/2,ind1:ind2);%%tfr(:,ind1:ind2);tfr(1:end/2,ind1:ind2);
% t1=t2(ind1:ind2);

% figure;
% %mesh(t1,f,tfr);
% pcolor(t1,f,tfr);
% shading interp;
% colorbar;
% xlabel('t/s');ylabel('f/Hz');title('振动信号的SPWV时频分布图')

% specnum=length(sig);
% 
% [tfr,rtfr,f] = tfrstft(sig,1:length(sig),length(sig),hamming(125),hamming(125));
% f=f*fs;
% figure;
% % mesh(t1,f,tfr);
% % pcolor(t2,f,abs(tfr));
% pcolor(t2,abs(f),abs(tfr));
% shading interp;
% colorbar;


clear ind1 ind2 fs specnum  tfr_stft; 
clear winlen  h vdata1 y f0;pack;
clear fp vdata1
%%% 出口：原来的fs_origin t vdata 及chishu参数
%%% 出口：变采后的时间fs,t2,sig  fp
%%% 出口：只是多加了两个截断后时间对齐标志ind1　ind2
%%
%*************%４）vitebi搜索频率路径****************************************
%%%入口参数：时间轴向量t1 频率轴向量f 时频分布tfr  起始频率startfreq 
%%%入口参数：原始振动信号的vdata,fs_origin ,t及齿数参数chishu 
% %%%%%%%划分频率轴格
df=(f(2)-f(1));meiMge=floor(1/df); 
[tfr2,~,ind3]=pinglvfengduan(tfr,f,meiMge);
clear tfr;tfr=tfr2;
clear tfr2 
%%%%%%%划分频率轴格
c=5;sigma=3;efxia=2;
[route,val_route]=viterbi3(abs(tfr),c,sigma);
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
f1=route(10,:);
%%%%对计算出来的估计转速插值回原来的fs
fs=fs_origin;
f1=interp1(t1,f1,t1(1):1/fs:t1(end));%%%%%
% f2=interp1(t1,f2,t1(1):1/fs:t1(end));
f1=f1/2;
t1=t1(1):1/fs:t1(end);%%%%%%
a=polyfit(t1,f1,3);
f2=polyval(a,t1,3);
figure;
hold on;
plot(t1,f1,'k');
plot(t1,f2,'r');
xlabel('t/s');ylabel('amplitude');title('估计转速');
hold off;
%%%%对计算出来的估计转速插值回原来的fs
%%
[~,ind1]=min(abs(t-jieduan1));[~,ind2]=min(abs(t-jieduan2)); 
vdata1=vdata(ind1:ind2);
% %%%%%%%%%%%%陷波滤波
vdata1=vdata1-mean(vdata1);
% vdata1=trapper(vdata1,50,fs,ceil(fs/50));
%%阶次分析
tiepian=1;
pluse=freq2pluse(f2,fs,tiepian);%%调入转频
figure(12)
hold on;
plot(t1,pluse);
ylim([0 1.3]);
% plot(pluse,'k')
xlabel('t/s');ylabel('amplitude');title('估计脉冲');
hold off;
order_max=40;
numpr=1;
threshold=abs(1)/3;
%%包络阶次分析
[order,ayy,array_angle_amp]=fuction_baoluopu(vdata1,pluse,fs,order_max,1,threshold);
figure;
plot(order,ayy);
title('包络阶次谱图'),
xlabel('order'),
ylabel('amplitude');
xlim([0.2 20]); 
%%阶次分析
[freq,ayy]=FUCTION_COTA2(vdata1,pluse,fs,order_max,1,threshold);
figure;
plot(freq,ayy)
title('阶次谱图'),
xlabel('order'),
ylabel('amplitude');
xlim([0.2 20]); 




