function[ex,ey,ra,dec,pa,tt]=read_dooley(froot)

fid=fopen([froot '_parallactic'],'r');
pa=fread(fid,inf,'double');
fclose(fid);

n=numel(pa);

fid=fopen([froot '_ra'],'r');
ra=fread(fid,inf,'double');
fclose(fid);
ndet=numel(ra)/n;
assert(ndet==round(ndet));

ra=reshape(ra,[ndet n])'*pi/180;

fid=fopen([froot '_dec'],'r');
dec=fread(fid,inf,'double');
fclose(fid);
dec=reshape(dec,[ndet n])'*pi/180;

fid=fopen([froot '_Ey2'],'r');
ey=fread(fid,inf,'double');
fclose(fid);
ey=reshape(ey,[ndet n])';

fid=fopen([froot '_Ex2'],'r');
ex=fread(fid,inf,'double');
fclose(fid);
ex=reshape(ex,[ndet n])';

try
  fid=fopen([froot '_time'],'r');
  tt=fread(fid,inf,'double');
  fclose(fid);
catch
  tt=(1:numel(pa))';
end


