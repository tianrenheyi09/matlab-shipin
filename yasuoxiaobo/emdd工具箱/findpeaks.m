function n = findpeaks(x)
% Find peaks. 找极大值点，返回对应极大值点的坐标
n    = find(diff(diff(x) > 0) < 0); % 一阶导大于零增函数，二阶导小于零凹函数
u    = find(x(n+1) > x(n));
n(u) = n(u)+1;                      % 加1才真正对应极大值点
