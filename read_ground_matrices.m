function[dat]=read_ground_matrices(fname)
fid=fopen(fname,'r');
dat.alt=fread(fid,1,'float64');
dat.az=fread(fid,1,'float64');
dat.ndata=fread(fid,1,'int');
dat.ndet=fread(fid,1,'int');
dat.rr=fread(fid,dat.ndet,'int');
dat.cc=fread(fid,dat.ndet,'int');
dat.match_list=fread(fid,dat.ndet,'int');
dat.std2=fread(fid,dat.ndet,'float64');
dat.std=fread(fid,dat.ndet,'float64');
dat.npoly=fread(fid,1,'int')+1; %+1 since legendre_mat returns 50th order poly, which is 51 terms
dat.daz=fread(fid,1,'float64');
dat.naz=fread(fid,1,'int');
dat.azpix=fread(fid,dat.naz,'int');
dat.rhs=fread(fid,[dat.npoly+dat.naz dat.ndet],'float64');
dat.mat=fread(fid,[dat.npoly+dat.naz dat.npoly+dat.naz],'float64');



cur_pos=ftell(fid);
info=stat(fid);
assert(cur_pos==info.size);
%disp(['read ' num2str(cur_pos) ' bytes out of ' num2str(info.size)]);
fclose(fid);

