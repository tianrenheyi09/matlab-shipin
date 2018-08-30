function [freq,ayy]=FUCTION_COTA3(vdata,pdata,Fs,order_max,Numpr,threshold)
%the main program of order tracking
% edit by Du zhao hui, 2011.12.12

% input parameters:
%vdata 输入的振动数据，当数据长度小于6e5时，参与计算的数据长度为整个数据长度；当数据长度大于6e5时，截取6e5数据分析
%pdata 输入的脉冲数据，数据为vdata的同步采样数据
%fs  采样频率
%order_max 需要分析的最大阶次（2<=order_max）
%Numpr     每转脉冲数目(轴上的键槽数目
%threshold  脉冲的阀值设定—脉冲最大值的2/3为阀值
%output parameters:
% 阶次谱图
% ---------------------------------------------------------------------
%                     Order Tracking MAIN PROCEDURE
% ----------------------------------------------------------------------
%取60000个数据点长度数据
numpr=1;
meiMpluse=Numpr;
if length(vdata)<=6e5
    array_time_amp=vdata;
    pluse=pdata;
else
    array_time_amp=vdata(1:6e5);
    pluse=pdata(1:6e5);
end
N=length(array_time_amp);       
t=(0:N-1)/Fs;
Num_pluse1=1;
for temp2=1:length(pluse)-1;
    if abs(pluse(temp2))<threshold&&abs(pluse(temp2+1))>=threshold
        Num_pluse1=[Num_pluse1,temp2+1];%计算脉冲数目
    end
end
a=[Num_pluse1(2):Num_pluse1(end)]
array_time_amp=array_time_amp(Num_pluse1(2):Num_pluse1(end));%脉冲点对应的幅值
Num_pluse1=Num_pluse1(2:end);
Num_pluse1=Num_pluse1-Num_pluse1(1)+1;
t=(0:length(array_time_amp)-1)/Fs;



if length(Num_pluse1)<2,return,end;
t_pluse=(Num_pluse1-1)/Fs;%脉冲点对应的时刻
 if mod(length(t_pluse),numpr)==0                         
     numpluse=floor((length(t_pluse)-1)/numpr);     
 else
     numpluse=floor(length(t_pluse)/numpr);
 end
%  if numpr>=2                                                          
%      whole_pluse=[];whole_pluse=[whole_pluse,t_pluse(1)];
%      for k=1:numpluse
%          imev=t_pluse(1+k*numpr);
%          whole_pluse=[whole_pluse,imev];%一圈只算一个脉冲时间点
%      end
%       t_pluse=whole_pluse;
%  else
      whole_pluse=t_pluse;
%  end

 len_whole_pluse=length(whole_pluse);
 lv=meiMpluse;
 num_len=ceil(len_whole_pluse/lv);
 whole_pluse2=[];
 for pp=1:num_len,
     whole_pluse2=[whole_pluse2,whole_pluse(1+(pp-1)*lv)];%一圈算同一个脉冲
 end
 whole_pluse=whole_pluse2;
i=1;
while i<=length(whole_pluse)-1;
    ft=1/(whole_pluse(i+1)-whole_pluse(i));%计算瞬时频率
     fomax=ft *order_max;%最大阶次对应的频率
      prodata=array_time_amp(Num_pluse1(1+(i-1)*numpr):(Num_pluse1(1+i*numpr)-1));
      prolen=length(prodata);fomn=ceil(prolen*fomax/Fs);
      freqd=fft(prodata);           
      if mod(prolen,2)==0         
          maxspectrumline=ceil(prolen/2)+1;      
      else
          maxspectrumline=ceil(prolen/2);         
      end
      if fomn<maxspectrumline;
          if mod(prolen-1,2)==0;          
              ii=maxspectrumline;
              while ii>fomn;
                  freqd(ii)=0;
                  freqd(mod(prolen-ii+1,prolen)+1)=0;
                  ii=ii-1;
               end
          else                                  
              freqd(maxspectrumline)=0;     
              ii=ceil(prolen/2);
              while ii>fomn;
                  freqd(ii)=0;
                  freqd(mod(prolen-ii+1,prolen)+1)=0;
                  ii=ii-1;
              end
          end;
      end
      Imed=ifft(freqd);
      array_time_amp(Num_pluse1(1+(i-1)*numpr):(Num_pluse1(1+i*numpr)-1))=Imed(1:(Num_pluse1(1+i*numpr)-Num_pluse1(1+(i-1)*numpr)));
      i=i+1;
end

delt_thet=pi/order_max;
t_angle=[];
for temp3=3:length(t_pluse);
    b=inv([1,t_pluse(temp3-2),t_pluse(temp3-2)^2;1,t_pluse(temp3-1),t_pluse(temp3-1)^2;1,t_pluse(temp3),t_pluse(temp3)^2])*[0,2*pi,4*pi]';
    if temp3==3;                                              
        k=0;
        while k<1.5*2*pi/delt_thet;
            if b(3)~=0;
                tt=(sqrt(4*b(3)*(k*delt_thet-b(1))+b(2)^2)-b(2))/(2*(b(3)+eps));
                t_angle=[t_angle,tt];
            else
                tt=(k*delt_thet-b(1))/b(2);
                t_angle=[t_angle,tt];
            end
            k=k+1;
        end
    else                                                             
        k=pi/delt_thet;
        while k>=pi/delt_thet && k<1.5*2*pi/delt_thet;
            if b(3)~=0;
                tt=(sqrt(4*b(3)*(k*delt_thet-b(1))+b(2)^2)-b(2))/(2*(b(3)+eps));
                t_angle=[t_angle,tt];
            else
                tt=(k*delt_thet-b(1))/b(2);
                t_angle=[t_angle,tt];
            end
            k=k+1;
        end
    end
end
k=3*pi/delt_thet;                                        
while k<4*pi/delt_thet;
    if b(3)~=0;
        tt=(sqrt(4*b(3)*(k*delt_thet-b(1))+b(2)^2)-b(2))/(2*(b(3)+eps));
        t_angle=[t_angle,tt];
    else
        tt=(k*delt_thet-b(1))/b(2);
        t_angle=[t_angle,tt];
    end
    k=k+1;
end
array_angle=[0:length(t_angle)-1].*delt_thet;
array_angle_amp=interp1(t,array_time_amp,t_angle,'spline');
array_angle_amp=array_angle_amp-mean(array_angle_amp);  

%时域平均
lv=10;
vec_num=floor(length(array_angle_amp)/(2*order_max*lv));
temp=zeros(1,2*order_max*lv);
for pp=1:vec_num,
   temp=temp+array_angle_amp((pp-1)*2*order_max*lv+1:(pp-1)*2*order_max*lv+2*order_max*lv);
end
array_angle_amp=temp;
array_angle_amp=array_angle_amp/vec_num;
array_angle=[0:length(array_angle_amp)-1].*delt_thet;


%%
angle_dom_ffty=abs(fft(array_angle_amp))*2/length(array_angle_amp);
 delt_order=2*pi/(length(angle_dom_ffty)*delt_thet);
angle_dom_fx=(0:length(angle_dom_ffty)-1)*delt_order;
%plot(angle_dom_fx(1:length(angle_dom_fx)/2),angle_dom_ffty(1:length(angle_dom_fx)/2)),
freq=angle_dom_fx(1:length(angle_dom_fx)/2);
ayy=angle_dom_ffty(1:length(angle_dom_fx)/2);

end
