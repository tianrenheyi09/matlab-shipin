clear all;
close all;
vib=textread('E:\bishe实验数据\25600数据\REC3939_ch3.txt','','headerlines',16);
pul=textread('E:\bishe实验数据\25600数据\REC3939_ch2.txt','','headerlines',16);
pluse1=vib(:,1);
vdata=vib(:,2);
% pul=pul(:,2);pul=-pul;
threshold=abs(1)/3;
fs=25600;
figure;
plot(pluse1,vdata)
xlabel('t/s');ylabel('amplitude');title('原始信号');
%%%取降速信号处理：
vdata=vdata(fs*7:fs*9);
% pul=pul(fs*7:fs*9);
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
%% 小波降噪
lev=5;
sig=wden(sig,'minimaxi','s','mln',lev,'sym5');%%陀螺方位角速度去噪
figure;
plot(t,sig);


% % 峭度滤波
% nlevel= 8;
% opt1=1;opt=2;
% figure()
% sig= Fast_Kurtogram2(sig,nlevel,fs,opt1,opt);
%% 低通滤波 
% fp=[5000,6100];fst=[4900,6200];%滤波器的参数设置


fp=[3800,4600];fst=[3600,4800];%滤波器的参数设置
% sig=sig-mean(sig);
% sig=Noisedaitong(sig,fs,fp,fst);
% figure;
% plot(t,sig);
% xlabel('t/s');ylabel('amplitude')
% title('滤波后信号')


% 
% figure;
% plot(t,sig);
% xlabel('t/s');ylabel('amplitude')
% title('滤波后信号')
% sig=sig-mean(sig); 
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
title('降频采样信号');
xlabel('t/s');ylabel('幅值');
fs=fs2; 
clear 　chongcaiyanlv  pp fst;clear chongcaiyanlv fs2;pack;

% x=sig1;
% CWTopt=struct('gamma',eps,'type','morlet','mu',6,'s',5,'om',0,'nv',64,'freqscale','linear'); %小波变换操作数
% 
% % CWT Synchrosqueezing transform
% [Tx, f1, Wx, as, Cw] = synsq_cwt_fw(t2, x, CWTopt.nv, CWTopt);
% xNew = synsq_cwt_iw(Tx, f1, CWTopt).';
% 
% figure;
% imagesc(t2, f1, abs(Tx));axis xy % CWT-SST时频图
% set(gca,'YDir','normal')
% colorbar;
% xlabel('时间 t/s');
% ylabel('频率 f/Hz');
% title('同步压缩时频图');

% figure;
% mesh(t2,f1,abs(Tx));

%% 连续小波变换时频图
s=sig1;
wavename='cmor3-3';
totalscal=2048;
Fc=centfrq(wavename); % 小波的中心频率
c=2*Fc*totalscal;
scals=c./(1:totalscal);
fff=scal2frq(scals,wavename,1/fs); % 将尺度转换为频率
coefs=cwt(s,scals,wavename); % 求连续小波系数
figure;
% f=(0:length(sig1)/2-1)*fs/length(sig1);
f=0.5*(0:length(sig1)-1)*fs/length(sig1);
imagesc(t2,fff,abs(coefs));
set(gca,'YDir','normal')
colorbar;
xlabel('时间 t/s');
ylabel('频率 f/Hz');
xlabel('时间 t/s');
ylabel('频率 f/Hz');
title('小波时频图');

% %% 局部最大值
% fw=fff;
% coes_y=coefs;
% N=length(s);
% a=zeros(1,N);
% b=zeros(1,N);
% b_p=zeros(1,N);
% [a(1),b(1)]=max(abs(coes_y(1:end,1)));
% temp1=30;%%开始搜索的频率
% b_p(1)=b(1)+temp1;
% for i=2:N
%     [a(i),b(i)]=max(abs(coes_y(temp1+b(i-1)-5:temp1+b(i-1)+5,i)));
%     temp1=b_p(i-1)-5-1;
%     b_p(i)=b(i)+temp1;
% %     coes_y(b_p(i)-5:b_p(i)+5,i)=0;
% end
% figure(1);
% 
% f1=fw(b_p);
% 
% % plot(t2,f1);
% p=polyfit(t2,f1,1);
% y0=polyval(p,t2);
% plot(t2,f1,t2,y0,'k');
% 
% 
% 
% fs=fs_origin;
% t1=t2;
% 
% f1=interp1(t1,f1,t1(1):1/fs:t1(end));%%%%%
% % f2=interp1(t1,f2,t1(1):1/fs:t1(end));
% f1=f1/2;
% t1=t1(1):1/fs:t1(end);%%%%%%
% 
% [ttt,frequency2]=shijispeed(pul,threshold,fs);
% 
% a=polyfit(t1,f1,1);
% f2=polyval(a,t1,1);
% figure;
% hold on;
% plot(t1,f1,'k');
% plot(t1,f2,'r');
% plot(ttt,frequency2);
% xlabel('t/s');ylabel('amplitude');title('估计转速');
% hold off;





%% viterbi搜索

tfr=abs(coefs);
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



route=route*(fs/2/totalscal);

%% 奇异点校正
[row,col]=size(route);
for i=1:row
    for j=2:col
        if route(i,j)-route(i,j-1)>=15
            route(i,j)=route(i,j-1);
        end
    end
end

f1=route(1,:);
figure;
plot(t2,f1);

fs=fs_origin;
t1=t2;

f1=interp1(t1,f1,t1(1):1/fs:t1(end));%%%%%
% f2=interp1(t1,f2,t1(1):1/fs:t1(end));
f1=f1/2;
t1=t1(1):1/fs:t1(end);%%%%%%

[ttt,frequency2]=shijispeed(pul,threshold,fs);

a=polyfit(t1,f1,1);
f2=polyval(a,t1,1);
figure;
hold on;
plot(t1,f1,'k');
plot(t1,f2,'r');
plot(ttt,frequency2);
xlabel('t/s');ylabel('amplitude');title('估计转速');
hold off;



%  for p=1:size(route,1),   
% route(p,:)=f(route(p,:));
% end
% f1=route(12,:);
% f1=route(16,:);3942 3
f1=route(7,:);%5%1-3秒的话是15，2-3的话是7

%%%%对计算出来的估计转速插值回原来的fs

fs=fs_origin;
f1=interp1(t1,f1,t1(1):1/fs:t1(end));%%%%%
% f2=interp1(t1,f2,t1(1):1/fs:t1(end));
f1=f1/2;
t1=t1(1):1/fs:t1(end);%%%%%%

[ttt,frequency2]=shijispeed(pul,threshold,fs);

a=polyfit(t1,f1,3);
f2=polyval(a,t1,3);
figure;
hold on;
plot(t1,f1,'k');
plot(t1,f2,'r');
plot(ttt,frequency2);
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
xlim([0 20]); 

%%阶次分析
[f,ayy,t_angle,array_time_amp]=FUCTION_COTA2(vdata1,pluse,fs,order_max,1,threshold);
figure;
plot(f,ayy)
title('阶次谱图'),
xlabel('order'),
ylabel('amplitude');
xlim([0 20]); 




