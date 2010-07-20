function[ra,dec]=read_pointing(fname)
fid=fopen(fname,'r');
assert(fid~=-1);
ndet=fread(fid,1,'int');
ndata=fread(fid,1,'int');
dat=fread(fid,[ndata*2 ndet],'float=>float');
ra=dat(1:ndata,:);
dec=dat(ndata+1:end,:);
