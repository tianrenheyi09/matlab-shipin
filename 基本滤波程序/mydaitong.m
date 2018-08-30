function y1= mydaitong(y,fp1,fp2,fst1,fst2,fs)
% fp1=250;fp2=320;
% fst1=200;fst2=350;
ws1=fst1*pi*2/fs;wp1=2*fp1*pi/fs;
ws2=fst2*pi*2/fs;wp2=2*fp2*pi/fs;
trw=min((wp1-ws1),(ws2-wp2));
M=ceil(6.2*pi/trw);
M=M+mod(M+1,2);
wc1=(ws1+wp1)/2;wc2=(wp2+ws2)/2;
fc1=wc1/pi;fc2=wc2/pi;
h1=fir1(M-1,[fc1,fc2],hanning(M)');

% figure;
% freqz(h1,1,512,fs);

y1=filter(h1,1,y);
% figure;
% plot(t,y1);
% n=length(y1);
% f=(0:n/2-1)*fs/n;
% s=abs(fft(y1-mean(y1)));
% figure;
% plot(f,s(1:n/2));