function y=Noisedaitong(sig,fs,fp,fst)
wp=fp*2/fs;
ws=fst*2/fs;
rp=0.1;
as=1;
b=design(fdesign.bandpass(fst(1),fp(1),fp(2),fst(2),as,rp,as,fs),'butter');
y=filter(b,sig);
end