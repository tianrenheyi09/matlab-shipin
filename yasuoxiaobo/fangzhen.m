
clear all;clc;
close all;
fs=1024;
t=0:1/fs:(5*1024-1)/fs;

% f1=5*t+2;f2=20+50*t;

% plot(t,f1,'r');hold on;plot(t,f2,'b');
sig=cos(2*pi*(2.5*t.^2+2*t))+sin(2*pi*(25*t.^2+20*t)); %另别的sig
%sig=cos(2*pi*(10*t.^2+2*t))+cos(2*pi*(20*t.^2+4*t)); %另别的sig
%sig=cos(2*pi*(120-100/3*t.^3))+cos(2*pi*(60-50/3*t.^3));
vdata=sig;
figure();plot(t,sig);

t=(0:length(sig)-1)/fs;
t2=t;

% figure;
% plot(t,vdata,'k');
% xlabel('t/s');ylabel('amplitude/g');title('变转速过程中的时域波形')


% nlevel= 2;
% opt1=2;opt2=1;
% figure()
% sig= Fast_Kurtogram2(sig,nlevel,fs,opt1,opt2);


%延拓
% sig=sig';
% [sig,startnode,xlength,tnode]=duichengyantuo(sig,fs);

%低通滤波
% chishu=1;
% fp=chishu*50;
% fst=chishu*50+10;
% sig=sig-mean(sig);
% sig=Noiseditong(sig,fs,fp,fst);
sig=sig-mean(sig);
% sig=sig(startnode:startnode+xlength-1);
% t2=t(tnode:tnode+xlength-1);
% sig=sig';
% figure;plot(t2,sig);
% title('滤波后的信号');xlabel('t/s');ylabel('幅值');

%降低采频
% fs2=512*2;chongcaiyanlv=floor(fs/fs2);
% fs_origin=fs;
% pp=1:(floor((length(sig)-1)/chongcaiyanlv)+1);
% pp=(pp-1)*chongcaiyanlv+1;
% sig=sig(pp);
% t2=t2(pp);
% figure;
% subplot(211);plot(t,vdata,'k');
% xlabel('t/s');ylabel('amplitude/g');title('变转速过程中的时域波形')
% subplot(212);plot(t2,sig,'k');
% xlabel('t/s');ylabel('amplitude/g');title('降频为1024Hz的采样信号')
% fs=fs2;
% clear fs2　chongcaiyanlv  pp fst;clear chongcaiyanlv fs2;pack;

%WVD去端点效应再截断两边各1s
jieduan1=1;jieduan2=4;
% sig=sig';
[~,ind1]=min(abs(t2-jieduan1));[~,ind2]=min(abs(t2-jieduan2));
sig=sig-mean(sig);
sig=sig';sig=sig-mean(sig);
specnum=length(sig);
sig=hilbert(sig);sig=abs(sig);sig=sig-mean(sig);
sig=hilbert(sig);
% g=tftb_window(length(sig)/10+1,'hamming'); h=tftb_window(length(sig)/10+1,'hamming');
% [tfr,~,f] = tfrspwv(sig,1:length(sig),specnum,g,h);
[tfr,~,f] = tfrspwv(sig,1:length(sig),specnum,hamming(251),hamming(251));
f=f*fs;
tfr=tfr(:,ind1:ind2);
t1=t2(ind1:ind2);
figure();
% mesh(t1,f,tfr);
pcolor(t1,f,tfr);
shading interp;
colorbar;
xlabel('t/s');ylabel('f/Hz');title('振动信号的SPWV时频分布图')

clear ind1 ind2 fs sig specnum t2 tfr_stft;
clear winlen  h vdata1 y f0;pack;
clear fp vdata1

%划分频率轴格
df=(f(2)-f(1));meiMge=floor(1/df);
[tfr2,~,ind2]=pinglvfengduan(tfr,f,meiMge);
clear tfr;tfr=tfr2;clear tfr2

%viterbi3算法
c=5;sigma=3;efxia=2;
[route,val_route]=viterbi3(abs(tfr),c,sigma,efxia);

%提取转频和二倍转频
for p=1:size(route,1),
    temp=route(p,:);
    temp2=zeros(1,length(temp));
    for time_node=1:length(temp),
        temp2(time_node)=(temp(time_node)-1)*meiMge+ind2(temp(time_node),time_node)-1;
    end
    route(p,:)=temp2;
end
clear temp temp2 time_node p;
for p=1:size(route,1),
    route(p,:)=f(route(p,:));
end
f1=route(1,:);
% f2=route(6,:);

%对计算出来的估计转速插值回原来的fs
fs=1024;
f1=interp1(t1,f1,t1(1):1/fs:t1(end));
f1=f1/10;
% f2=interp1(t1,f2,t1(1):1/fs:t1(end));
% f2=f2/3;
t1=t1(1):1/fs:t1(end);
figure();
% hold on;
plot(t1,f1,'k-');
xlabel('t/s');ylabel('f/Hz');title('Viterbi算法提取的转速')
% plot(t1,f2,'b-.');
% hold off;


[~,ind1]=min(abs(t-1));[~,ind2]=min(abs(t-4));
vdata1=vdata(ind1:ind2);

%陷波滤波
vdata1=vdata1-mean(vdata1);
vdata1=trapper(vdata1,50,fs,ceil(fs/50));

vdata1=vdata1-mean(vdata1);
% vdata1=hilbert(vdata1);
% vdata1=abs(vdata1);
% vdata1=vdata1-mean(vdata1);

tiepian=1;
pluse=freq2pluse(f1,fs,tiepian);
figure; 
plot(t1,pluse,'k');ylim([0 1.3]);xlabel('t/s');ylabel('amplitude');title('键相信号的估计')
order_max=40;numpr=1;threshold=abs(1)*1/3;
[f,ayy]=FUCTION_COTA(vdata1,pluse,fs,order_max,1,threshold);
figure();
plot(f,ayy,'k');xlabel('t/s');ylabel('amplitude');title('基于故障特征频率的包络阶次分析图')
%}