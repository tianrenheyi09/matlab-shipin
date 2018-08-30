function [A,f,tt] = hhspectrum(x,t,l,aff)
error(nargchk(1,4,nargin));
if nargin < 2
t=1:size(x,2);
end

if nargin < 3
l=1;
end

if nargin < 4
aff = 0;
end

if min(size(x)) == 1
if size(x,2) == 1
x = x';
if nargin < 2
t = 1:size(x,2);
end
end
Nmodes = 1;
else
Nmodes = size(x,1);
end

lt=length(t);
tt=t((l+1):(lt-l));

for i=1:Nmodes
an(i,:)=hilbert(x(i,:)')';
f(i,:)=instfreq(an(i,:)',tt,l)';
A=abs(an(:,l+1:end-l));

if aff
disprog(i,Nmodes,max(Nmodes,100))
end

end
