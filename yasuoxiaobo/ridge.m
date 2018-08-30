function out=ridge ( c ) 
  %c是要提取脊线的时频矩阵,c元素为非负
[B,A]=size(c);
N=floor(A*B/4);
%生成N个Climber
temp=zeros(B,A);
for i=1:1:A*B
   if mod(i,4)==0
      temp(i)=1;
   end
end
T=max(max(c))-min(min(c));%系统初始化温度
Tt=T; %系统当前温度
t=2; %系统当前时间

while Tt>=T/1000 %对时间t做循环
   for i=4:4:A*B %对每个climber做移动
      if mod(i,B)==0
         heng=mod(i,B)+B;
      else
         heng=mod(i,B);
      end %计算climber的横坐标
     
      zong=ceil(i/B) %计算climber的纵坐标
      p=sign(2*rand-1); %横坐标以0.5的等概率分别向左和右移动
      if heng==1
         p=1;
      elseif heng==A
         p=-1;
      end %排除边界条件
      heng_new=heng+p;

      %纵坐标按规则移动
      p=sign(2*rand-1);
      if zong==1
         p=1;
      elseif zong==A
         p=-1;
      end %排除边界条件
      zong_new=zong+p;
      if c(heng_new,zong_new)>c(heng_new,zong)
         temp(heng,zong)=0;
         temp(heng_new,zong_new)=1;
      else
         pt=exp((c(heng_new,zong_new)-c(heng_new,zong))/Tt);
         if(rand<=pt)
             zong_new=zong+p;
             temp(heng,zong)=0;
             temp(heng_new,zong_new)=1;
         else
             zong_new=zong;
             temp(heng,zong)=0;
             temp(heng_new,zong_new)=1;
          end 
      end
   end
   Tt=T/log2(t);
   t=t+1;
end
