function [ttt,frequency2]=shijispeed(pul,threshold,fs)
% pul=textread('G:\code_for1\REC3950_ch2.txt','','headerlines',16);
% pul=pul(:,2);
% pul=-pul;
% fs=25600;
% threshold=abs(1)/3;
% pul=pul(fs*11:fs*16);


N=length(pul);       
t=(0:N-1)/fs;
time=t;
Num_pluse1=1;

k_p=find(pul<0);
pul(k_p)=0;
pul=pul/max(pul);
Num_pluse1=1;
for temp2=1:length(pul)-1;%去除负脉冲
    if abs(pul(temp2))<threshold&&abs(pul(temp2+1))>=threshold
        Num_pluse1=[Num_pluse1,temp2+1];
    end
end
k=[];
for i=2:length(Num_pluse1)
    if (Num_pluse1(i)-Num_pluse1(i-1))<256
        k=[k,i];
    end
end
m=0;
for i=1:length(k)
    Num_pluse1(k(i)-m)=[];
    m=m+1;
end
pul(Num_pluse1)=1;
% figure(3); 
% % subplot(211);
% stem(t(Num_pluse1),pul(Num_pluse1));ylim([0 1.3]);
% xlabel('t/s');ylabel('amplitude');title('处理后的脉冲信号');


% figure(4); 
% n=length(vib);
% f_vdata=abs(fft(vib))*2/n;f=(0:n-1)*fs/n;
% plot(f,f_vdata);
% xlabel('frquency/Hz');ylabel('amplitude/g');title('频谱图');
% xlim([0 fs/2]);
% %ylim([0 0.12]);

frequency2=zeros(1,length(Num_pluse1)-2);
ttt=frequency2;
for i=3:length(Num_pluse1)
    frequency2(i-1)=1/(time(Num_pluse1(i))-time(Num_pluse1(i-1)));
    ttt(i-1)=time(Num_pluse1(i));
end
ttt(1)=[];frequency2(1)=[];
% tttt=tt;
figure;
plot(ttt,frequency2);

xlabel('t/s');ylabel('Fre/Hz');title('估计转频');
end




