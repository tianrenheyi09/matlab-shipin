function pluse=freq2pluse(f,fs,tiepian)
%[R]=romberg(f,a,b,fs,tol)

%%%1)设置参数
pluse=zeros(1,length(f));
tol=1e-6;
round=tiepian;%round标记是几转,tiepian是模拟轴上贴几个反光片
%%%1)设置参数

%%%2）减少计算量，先用梯形法粗略积分，再在粗略点附近用romberg积分法
R1=0;gujishibiao=[];r1=[];
round=1;
for num=1:length(f)-1,
    R1=R1+(f(num)+f(num+1))/fs/2;
    if R1>round,
        gujishibiao=[gujishibiao,num];
        round=round+1;
        r1=[r1,R1];
    end
end
R1=r1;
%%%2）减少计算量，先用梯形法粗略积分，再在粗略点附近用romberg积分法

%%%3）对估计时标处用romberg方法再积分，得到较精确的值R1
% R1=[];
% for p=1:length(gujishibiao),
%     a=0;b=gujishibiao(p)/fs;
%     temp=romberg(f,a,b,fs,tol);
%     R1=[R1,temp];
% end
%%%3）对估计时标处用romberg方法再积分，得到较精确的值R1

%%%4）对估计点处的值用romberg方法再算一遍，记录偏大偏小指标updown
round=1;updown=[];up=1;down=0;
for p=1:length(R1),
    if R1(p)>round,
        updown=[updown,up];
    else
        updown=[updown,down];
    end
    round=round+1;
end
%%%4）对估计点处的值用romberg方法再算一遍，记录偏大偏小指标updown

num=1;round=tiepian;ind=1;
up_biaozhi=1;down_biaozhi=-1;

%%%5）对R1结果修正，得到R和新的时标位置指针ind
ind=2; %键相标在时间序列中的指针，第一个是上升脉冲，标记以第1个数据定义初始键相
for p=1:length(updown),
    if updown(p)==1,
        %%%如是up可以表明计算误差使计算结果偏大，要用romberg方法反
        %%%反向搜索几步，以减少误差
        while up_biaozhi>0,
            num1=gujishibiao(p);
            num2=num1-1;            
            zhengshu_f=ceil(f(num1));  %gujishubiao向下处的整数值
            temp=(f(num1)+f(num2))/fs/2;
            R1(p)=R1(p)-temp;
            up_biaozhi=R1(p)-zhengshu_f;
            num1=num1-1;num2=num1-1;
        end
    else
        %%%如是down可以表明计算误差使计算结果偏小，要用romberg方法向前
        %%%正向搜索几步，以减少误差  
        while down_biaozhi<0,
            num1=gujishibiao(p);
            num2=num1+1;            
            zhengshu_f=floor(f(num1));  %gujishubiao向下处的整数值
            temp=(f(num1)+f(num2))/fs/2;
            R1(p)=R1(p)+temp;
            down_biaozhi=R1(p)-zhengshu_f;
            num1=num1+1;num2=num1+1;    
        end
    end
    ind=[ind,num1];
    up_biaozhi=1;down_biaozhi=-1;
end
%%%5）对R1结果修正，得到R和新的时标位置指针ind

% for num=1:length(f)-1,
%     temp_t=t(num);
%     a=0;b=temp_t;
%     [R1]=romberg(f,a,b,fs,tol);
%     b2=t(num+1);
%     [R2]=R1+(f(num)+f(num+1))/fs/2;
%     if  R1<=round&R2>=round,
%         ind=[ind,n];
%         round=round+1;
%     end
% end
pluse(ind)=1;
end