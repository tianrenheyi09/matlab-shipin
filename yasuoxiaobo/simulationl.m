clc;clear;close all
t1=0:1/40:7.5-1/40;
t2=7.5:1/40:25.6-1/40;
t=[t1,t2];
x1=sin(1.6*2*pi*t1)+sin(3*2*pi*t1+0.04*2*pi*t1.^2+4*sin(0.1*2*pi*t1));
x2=sin(1.6*2*pi*(t2-7.5)+0.055*2*pi*(t2-7.5).^2)+sin(3*2*pi*t2+0.04*2*pi*t2.^2+4*sin(0.1*2*pi*t2));
x=[x1,x2];
N=length(x);
y=awgn(x,11);
fs=40;fc=3;totalscal=200;
FreqBins=linspace(1/fs,6/fs,totalscal);
fw=FreqBins*fs;
scales=fc./FreqBins;

coes_y=cwt(y,scales,'cmor 4-3');
imagesc(t,scales,abs(coes_y));
figure(2)
imagesc(t,fw,abs(coes_y));
a=zeros(1,N);
b=zeros(1,N);
b_p=zeros(1,N);


[a(1),b(1)]=max(abs(coes_y(81:121,1)));
temp1=80;b_p(1)=b(1)+temp1;
for i=2:N
    [a(i),b(i)]=max(abs(coes_y(temp1+b(i-1)-10:temp1+b(i-1)+10,i)));
    temp1=b_p(i-1)-10-1;
    b_p(i)=b(i)+temp1;
end
figure(3)
plot(t,fw(b_p));
% hold on;

c=zeros(1,N);
d=zeros(1,N);
d_p=zeros(1,N);
[c(1),d(1)]=max(abs(coes_y(11:41,1)));
temp2=10;d_p(1)=d(1)+temp2;
for i=2:N
   
    [c(i),d(i)]=max(abs(coes_y(temp2+d(i-1)-10:temp2+d(i-1)+10,i)));
    
    temp2=d_p(i-1)-10-1;
    d_p(i)=d(i)+temp2;
    
end
figure(3)
hold on
plot(t,fw(d_p),'r')

F1=zeros(21,N);
e1=zeros(1,N);
p1=zeros(1,N);
pp1=zeros(1,N);

lamda=2;
F1(:,2)=-abs((coes_y(b_p(1),1)^2)+coes_y(b_p(2),2)^2)+lamda*(fw(b_p(2))-fw(b_p(1)))^2;
for j=3:N
for i=1:21
    F1(i,j)=F1(1,2)-abs(coes_y(b_p(j)+i-11,j)^2)+lamda*(fw(b_p(j)+i-11)-fw(b_p(j-1)+pp1(j-1)))^2;
end
[e1(j),p1(j)]=min(F1(:,j));
pp1(j)=p1(j)-11;
end
plot(t,fw(b_p+pp1),'g')

F2=zeros(21,N);
e2=zeros(1,N);
p2=zeros(1,N);
pp2=zeros(1,N);


F2(:,2)=-abs((coes_y(d_p(1),1)^2)+coes_y(d_p(2),2)^2)+lamda*(fw(d_p(2))-fw(d_p(1)))^2;
for j=3:N
for i=1:21
    F2(i,j)=-abs(coes_y(d_p(j)+i-11,j)^2)+lamda*(fw(d_p(j)+i-11)-fw(d_p(j-1)+pp2(j-1)))^2;
end
[e2(j),p2(j)]=min(F2(:,j));
pp2(j)=p2(j)-11;
end
plot(t,fw(d_p+pp2),'g')








