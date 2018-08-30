function [tfr2,f2,ind2]=pinglvfengduan(tfr,f,meiMge)
[row,col]=size(tfr);
max_start=floor((row-1)/meiMge);
tfr2=zeros(max_start,col);ind2=zeros(max_start,col);
p=1:max_start;
start=(p-1)*meiMge+1;
f2=[];
for p=1:max_start,
    temp_start=start(p);temp_end=start(p)+meiMge-1;
    temp_tfr=tfr(temp_start:temp_end,:);
    [temp_val,temp_ind]=max(temp_tfr);%返回最大和索引值搜索每一列的最大值然后返回第几行
    tfr2(p,:)=temp_val;%每一列最大值，一共col列
    ind2(p,:)=temp_ind;%ind2没一行对应10行中每一列最大值对应的第几行
    f2=[f2,f(temp_start)];%将f2化为1，11，21
end
f2=f2';
end