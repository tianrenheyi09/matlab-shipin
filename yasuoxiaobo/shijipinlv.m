function [ttt,frequency2]=shijipinlv(pul,fs,m,n)

% fs=12800;

threshold=abs(1)*1/3;

% pul=textread('E:\viterbi+阶次仿真\REC3939_ch2.txt','','headerlines',16); %3940为1寸轴承外圈信号



time=pul(:,1);
pul=pul(:,2);
pul=-pul;







time=time(fs*m:fs*n);%取升速信号

pul=pul(fs*m:fs*n);






N=length(pul);       
t=(0:N-1)/fs;
time=t;
% figure(2); 
% % subplot(211);
% plot(t,pul);xlabel('t/s');ylabel('amplitude');title('截取脉冲信号');

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
scatter(ttt,frequency2,'*');
hold on;
plot(ttt,frequency2);
hold off;
xlabel('t/s');ylabel('Fre/Hz');title('估计转频');
end




