clear;clc;close all
SamFreq=2048;
N=2*SamFreq;
fs=SamFreq;
t=0:1/fs:(N-1)/fs;

s_f1=(2.5*t.^2+10*t);s_f2=(12.5*t.^2+50*t);
f1=5*t+10;f2=50+25*t;

% s_f1=(5*t.^2+10*t);s_f2=(25*t.^2+50*t);
% f1=10*t+10;f2=50+50*t;

f1_zheng=f1;f2_zheng=f2;
figure(1);
plot(t,f1,'r');hold on;plot(t,f2,'b');
% sig=cos(2*pi*s_f1)+sin(2*pi*s_f2); %另别的sig
sig=cos(2*pi*s_f1)+sin(2*pi*s_f2); %另别的sig
%sig=cos(2*pi*(10*t.^2+2*t))+cos(2*pi*(20*t.^2+4*t)); %另别的sig
%sig=cos(2*pi*(120-100/3*t.^3))+cos(2*pi*(60-50/3*t.^3));

% sig=awgn(sig,-2);

 sig=sig';
 SNR=1;
 sig=sigmerge(sig,rand(N,1),SNR);%sigmerge要使sig和rand（N,1）保持一列，也就是sigmerge函数是N行一列
vdata=real(sig);
%vdata=sig;
figure(2)
plot(t,vdata);
xlabel('t/s');
ylabel('A/g');
set(gcf,'unit','centimeters','position',[10 8 13 9])


%%频率估计
sig=vdata;
t=(0:length(sig)-1)/fs;

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



%% spwvd
% [tfr,~,f] = tfrspwv(sig,1:length(sig),length(sig),hamming(125),hamming(125));
% [tfr,~,f] = tfrwv(sig);
% f=f(1:end/2)*fs;%%%%
% tfr=tfr(1:end/2,:);%%tfr(:,ind1:ind2);tfr(1:end/2,ind1:ind2);
% t1=t2;
% figure;
% %mesh(t1,f,tfr);
% pcolor(t1,f,tfr);
% shading interp;
% colorbar;
% xlabel('t/s');
% ylabel('f/Hz');
% set(gcf,'unit','centimeters','position',[10 8 13.5 9])
% title('WVD时频图')


% sig=hilbert(sig);
% sig=sig-mean(sig);
% % [tfr,~,f] = tfrspwv(sig,1:length(sig),specnum,g,h);
% [tfr,rtfr,f] = tfrstft(sig,1:length(sig),length(sig),hamming(65));
% tfr=tfr(1:end/2,:);
% f=f(1:end/2,:)*fs;
%  figure();
%  imagesc(t2, f, abs(tfr));
% set(gca,'YDir','normal')
% colorbar;
% xlabel('t/s');
% ylabel('f/Hz');
% title('短时傅里叶时频图')
% set(gcf,'unit','centimeters','position',[10 8 13 9])



% % % 连续小波
% s=sig;
% wavename='cmor3-3';
% % wavename='shan2-3';
% % wavename='mexh';
% totalscal=128;
% Fc=centfrq(wavename); % 小波的中心频率
% c=2*Fc*totalscal;
% scals=c./(1:totalscal);
% f=scal2frq(scals,wavename,1/fs); % 将尺度转换为频率
% coefs=cwt(s,scals,wavename); % 求连续小波系数morlet
% figure;
% f=(0:length(sig)/2-1)*fs/length(sig);
% 
% imagesc(t2,f,abs(coefs));
% set(gca,'YDir','normal')
% colorbar;
% xlabel('t/s');
% ylabel('f/Hz');
% title('小波时频图');
% set(gcf,'unit','centimeters','position',[10 8 13 9])


x=sig;
CWTopt=struct('gamma',eps,'type','hshannon','mu',6,'s',5,'om',0,'nv',128,'freqscale','linear'); %小波变换操作数

% CWT Synchrosqueezing transform
[Tx, f1, Wx, as, Cw] = synsq_cwt_fw(t2, x, CWTopt.nv, CWTopt);
xNew = synsq_cwt_iw(Tx, f1, CWTopt).';

figure;
imagesc(t2, f1, abs(Tx));axis xy % CWT-SST时频图
set(gca,'YDir','normal')
colorbar;
xlabel('t/s');
ylabel('f/Hz');
title('同步压缩时频图');
set(gcf,'unit','centimeters','position',[10 8 13 9])



f=f1;
tfr=abs(Tx);
df=(f(2)-f(1));meiMge=floor(1/df); 
[tfr2,~,ind2]=pinglvfengduan(tfr,f,meiMge);
clear tfr;tfr=tfr2;clear tfr2 
% %%%%%%%划分频率轴格

c=5;sigma=3;efxia=2;
[route,val_route]=viterbi3(abs(tfr),c,sigma,efxia);
for p=1:size(route,1),   
    temp=route(p,:);
    temp2=zeros(1,length(temp));
    for time_node=1:length(temp),
        temp2(time_node)=(temp(time_node)-1)*meiMge+ind2(temp(time_node),time_node)-1;
    end
    route(p,:)=temp2;
end
clear temp temp2 time_node p;
% route=route*(fs/2/totalscal);
% ff=route;%%搜索的频率，注意和stft与spwv搜的区别

 for p=1:size(route,1),   
route(p,:)=f(route(p,:));
 end
% route=route*(fs/2/);
f1=route(1,:);
f2=route(2,:);

jieduan1=0.5;jieduan2=1.5
[~,ind1]=min(abs(t2-jieduan1));[~,ind2]=min(abs(t2-jieduan2)); 
f1=f1(ind1:ind2);
f2=f2(ind1:ind2);























% %%
% jieduan1=1;jieduan2=ceil(max(t2)-1);
% %*************%２）计算WVD后要截取的信号指针****************************************
% %%% 入口：变采后的时间fs,t2,sig  fp
% %%% 入口：原来的fs_origin t vdata 及chishu参数
% %%% 目的：在信号两端留下１s时长以便清除WVD端点效应
% sig=sig';
% %%%%WVD去端点效应再截断两边各1s,实际上是第jieduan1+1:jieduan2-1(2~5s)的信号
% [~,ind1]=min(abs(t2-jieduan1));[~,ind2]=min(abs(t2-jieduan2));
% sig=sig-mean(sig);
% % %%%%%%%A)stft
%  sig=sig';sig=sig-mean(sig);
%  specnum=length(sig);
%  %sig=hilbert(sig);sig=abs(sig);sig=sig-mean(sig);
%  sig=hilbert(sig);
% g=tftb_window(length(sig)/10,'hamming'); h=tftb_window(length(sig)/10,'hamming'); 
% 
% if mod(length(g),2)==0,
%     g=tftb_window(length(sig)/10+1,'hamming');
% end
%    
% if mod(length(h),2)==0,
%    h=tftb_window(length(sig)/10+1,'hamming');
% end
% %% stft+spwv
% % [tfr,~,f] = tfrspwv(sig,1:length(sig),specnum,g,h);
% [tfr,rtfr,f] = tfrstft(sig,1:length(sig),specnum,g,h);
% %%%%%%%D)
% f=0.5*(0:length(sig)-1)*fs/length(sig)'
% % f=f*fs;
% % DD=(0.5*(0:length(f)-1)*fs/length(f))';
% tfr=tfr(:,ind1:ind2);
% 
% t1=t2(ind1:ind2);
% 
%  figure();
%  imagesc(t, f, abs(tfr));
% set(gca,'YDir','normal')
% colorbar;
% %mesh(t1,f,tfr);
% % pcolor(t1,abs(f),abs(tfr));
% % shading interp;
% % colorbar;
% xlabel('t/s');ylabel('f/Hz');title('振动信号的SPWV时频分布图')
% set(gcf,'unit','centimeters','position',[3 5 13.5 9])
% set(gcf,'color','white');
% ylim([0 200]);



%%%%对计算出来的估计转速插值回原来的fs
t1=t2;
t1=t1(ind1:ind2);
fs=fs_origin;
f1=interp1(t1,f1,t1(1):1/fs:t1(end));
f2=interp1(t1,f2,t1(1):1/fs:t1(end));
t1=t1(1):1/fs:t1(end);


figure(1);hold on;
plot(t1,f1,'b-','linewidth',2);
plot(t1,f2,'k','linewidth',2);hold off;

set(gcf,'unit','centimeters','position',[10 8 13 9])
xlabel('t/s'),
ylabel('f/hz');
title('对比图');

%%%%对计算出来的估计转速插值回原来的fs

%%%出口参数：时间轴向量t2 频率轴向量f   
%%%出口参数：原始振动信号的vdata,fs_origin ,t及齿数参数chishu 
%*************%４）vitebi搜索频率路径****************************************
%%

%%
[~,ind1]=min(abs(t-jieduan1));[~,ind2]=min(abs(t-jieduan2)); 
vdata1=vdata(ind1:ind2);
vdata1=vdata1-mean(vdata1);

% %%
% f1_zheng=f1_zheng(ind1:ind2);
% f2_zheng=f2_zheng(ind1:ind2);
% 
% per1=percent_err(f1,f1_zheng)
% per2=percent_err(f2,f2_zheng)

%%

%%
%f2=f2/32;
tiepian=1;
pluse=freq2pluse(f1,fs,tiepian);
figure;
plot(t1,pluse);
ylim([0,1.2]);

xlabel('t/s'),
ylabel('A');
title('估计的键相信号'),
set(gcf,'unit','centimeters','position',[10 8 13 9])

order_max=30;numpr=1;threshold=abs(1)*2/3;
%COTA(vdata,pdata,Fs,order_max,numpr,threshold)
figure;
[f,ayy]=FUCTION_COTA3(vdata1,pluse,fs,order_max,numpr,threshold);
plot(f,ayy);
title('阶次分析'),
xlabel('阶次'),
ylabel('幅值');
set(gcf,'unit','centimeters','position',[10 8 13 9])