function disp_hhs(im,t,inf)
% DISP_HHS(im,t,inf)
% displays in a new figure the spectrum contained in matrix "im"
% (amplitudes in log).
%
% inputs : - im : image matrix (e.g., output of "toimage")
% - t (optional) : time instants (e.g., output of "toimage") 
% - inf (optional) : -dynamic range in dB (wrt max)
% default : inf = -20
%
% utilisation : disp_hhs(im) ; disp_hhs(im,t) ; disp_hhs(im,inf) 
% disp_hhs(im,t,inf)

figure
colormap(bone)
colormap(1-colormap);

if nargin==1
inf=-20;
t = 1:size(im,2);

end

if nargin == 2
if length(t) == 1
inf = t;
t = 1:size(im,2);
else
inf = -20;
end
end

if inf >= 0
error('inf doit etre < 0')
end

M=max(max(im));

im = log10(im/M+1e-300);

inf=inf/10;


imagesc(t,fliplr((1:size(im,1))/(2*size(im,1))),im,[inf,0]);
set(gca,'YDir','normal')
xlabel(['time'])
ylabel(['normalized frequency'])
title('Hilbert-Huang spectrum')
