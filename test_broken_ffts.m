ndata=131000;nn=25;nn2=20
fftw_init_threaded
tod=allocate_tod_c();set_tod_dt_c(tod,1/200);set_tod_ndata_c(tod,ndata);
rows=1:nn;rows=repmat(rows,[nn2 1]);
cols=1:nn2;cols=repmat(cols,[nn 1])';
rows=reshape(rows,[nn*nn2 1]);
cols=reshape(cols,[nn*nn2 1]);set_tod_rowcol_c(tod,rows,cols);
allocate_tod_storage(tod);dat=get_tod_data(tod);dat=randn(size(dat));push_tod_data(dat,tod);clear dat
debutterworth_c(tod)
free_tod_storage(tod);