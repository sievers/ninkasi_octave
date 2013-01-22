function[value]=read_tod_data_new(tod)
allocate_tod_storage(tod);
fname=get_tod_name(tod);
[rr,cc]=get_tod_rowcol(tod);
dat=get_tod_data(tod);
ndet=length(rr);

myf=init_getdata_file(fname);

for j=1:ndet
  channame=sprintf('tesdatar%02dc%02d',rr(j),cc(j));
  vec=getdata_double_channel(myf,channame);
  dat(:,j)=vec;
end
push_tod_data(dat,tod);

close_getdata_file(myf);
