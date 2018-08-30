function y1= mygaotong(y,fp,fst,fs)
wp=2*fp*pi/fs;ws=2*fst*pi/fs;
deltaw=abs(ws-wp);
N=ceil(6.6*pi/deltaw);
N=N+mod(N,2);
wind=(hamming(N+1))';
Wn=(fp+fst)/fs;
% b=fir1(N,Wn,wind,'high');
b=fir1(N,Wn,'high');%% b=fir1(50,Wn,'high');
% figure;
% freqz(b,1,512,fs);


y1=filter(b,1,y);
% figure;
% plot(t,y1);

% n=length(y1);
% f=(0:n/2-1)*fs/n;
% s=abs(fft(y1-mean(y1)));
% figure;
% plot(f,s(1:n/2));
