function y=Noiseditong(sig,fs,fp,fst)
y=[];
%fp=1500;fst=2250;fs=15000; %输入设计指标
wp=2*fp/fs;                %求归一化数字通带截止频率
ws=2*fst/fs;                %求归一化数字阻带起始频率
% deltaw=ws-wp;               %求过渡带宽
% N0=ceil(6.6/deltaw);       %求窗口长度    
% N=N0+mod(N0+1,2);           %确保窗口长度N为奇数
% n=N-1;                      %求出滤波器的阶数n
n=100;
wn=(ws+wp)/2;                %求出滤波器的截止频率
b=fir1(n,wn);              %利用fir1函数求出滤波器的系数
% y=filter(b,1,sig);
% figure;
% freqz(b,1,512,fs);
% [h,w]=freqz(b,1,512,fs);
% magH=abs(h);
% phaH1=unwrap(angle(h));
% phaH=180*phaH1/pi;
% wnyq=w/pi;
% figure;
% subplot(2,1,1);
% plot(wnyq,magH);
% xlabel('频率f'); ylabel('幅值(dB)'); grid;
% subplot(2,1,2);
% plot(wnyq,phaH);
% xlabel('频率f'); ylabel('相位(degree)'); grid;	

%figure()
 y=filtfilt(b,1,sig);


end