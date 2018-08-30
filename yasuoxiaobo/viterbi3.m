function [route,val_route]=viterbi3(tfr,c1,sigma0,efxia)

[row,col]=size(tfr);

tfr(1:3,:)=-1000*ones(3,col);

[row,col]=size(tfr);
sigma=sigma0;

temp=zeros(row,1);
%   N1:每列取的搜索点数
%   N2:每个最大点搜索邻近的点数
N1=20;N2=20*4;
% pre_route=inf*ones(row,col);
% val_pre_route=inf*ones(row,col);

%%%%%%%记算代价函数f(x)=trf(x)
    for di_i=1:col,
        di_i_col=tfr(:,di_i);
        [val,ind]=sort(di_i_col,'descend');
        temp1=[0:row-1]';
        temp(ind)=temp1;
        tfr(:,di_i)=temp;
    end
%%%%%%%记算代价函数f(x)=tfr(x)
   
%%%route各列记录是网格中绝对网格点
route=inf*ones(N1,col);

[c_mins,ind_mins]=sort(tfr(:,1));
ind_mins=ind_mins(1:N1);%indmins是从小到大排列的数在原矩阵的第几行
route(:,1)=ind_mins;%route行数为n1

val_route=c_mins(1:N1);

%%%%%%%计算每点的后向最佳路径pre_route及val_pre_route
  for di_i=1:col-1,    %计算di_i列
%       [c_mins,ind_mins]=sort(tfr(:,di_i));
%       ind_mins=ind_mins(1:N1);
%       route(:,di_i)=ind_mins;
        ind_mins=route(:,di_i);
    for pointer=1:N1,%比较di_i列的第pointer个点
        %向后搜索
          min_pre_ind=max([ind_mins(pointer)-N2,1]);
          max_pre_ind=min([ind_mins(pointer)+N2,row]);
          keneng_ind=[min_pre_ind:max_pre_ind];
%        pre_pinglv=ind_mins(pointer);%保留前一个频率结点,pointer行，di_i列%%
          %%%%%%%%%%%%%%%计算代价函数g(x)
          gx=abs(keneng_ind-ind_mins(pointer));
          ind_sigma=find(gx<sigma|gx==sigma);
          ind_sigma2=find(gx>sigma);
          gx(ind_sigma)=0;
          gx(ind_sigma2)=c1*(gx(ind_sigma2)-sigma);
          %%%%%%%%%%%%%%%计算代价函数g(x)
          temp=gx'+tfr(keneng_ind,di_i+1);  
          [c0,ind0]=min(temp);  %ind0代表di_i-1列中keneng_ind（ind0）结点
                                %结点tfr(keneng_ind(ind0),di_i-1)
                                %与本结点tfr(ind_mins(pointer),di_i)最近
                                %不过可能不是前面搜出的较好的结点
          chang_ind0=find(temp==c0); %chang_ind0是前向代价函数中最小的一些
                                     %意义：tfr(keneng_ind(chang_ind0,di_i))
          foundind=keneng_ind(chang_ind0); %foundind是全局的
          [found_c,found_ind]=min(abs(foundind-ind_mins(pointer)));
          route(pointer,di_i+1)=foundind(found_ind); 
          val_route(pointer)=val_route(pointer)+c0;
%           xianzaipinglv=foundind(found_ind);%现在频率结点,pointert行，di_i+1列
%                     number=0;
%                     if xianzaipinglv==pre_pinglv,
%               number=number+1;
%               if number>20,
%                   sigma=sigma+0.1*(number-20);
%                   pre_pinglv=xianzaipinglv;
%               end
%           else
%               sigma=sigma0;
%               number=0;
%           end
    end
  end
end