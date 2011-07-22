function[timestreams]=get_bad_timestreams_mustang(tods,myopts)

freq=get_struct_mem(myopts,'pulse_freq',1.411);
nufrac=1/get_struct_mem(myopts,'pulse_oversample',2);
nnuvec=get_struct_mem(myopts,'pulse_nvec',1);
npoly=get_struct_mem(myopts,'npoly_bad',1);
for j=1:length(tods),
  [rows,cols]=get_tod_rowcol(tods(j));
  vecs_pt=get_pulse_tube_vecs(tods(j),freq,nufrac,nnuvec);
  %vecs_poly=legendre_mat(get_tod_ndata(tods(j)),npoly);
  %vecs_poly=vecs_poly(:,2:end);
  vv=(1:get_tod_ndata(tods(j)))';
  vv=vv-mean(vv);
  vv=vv/max(vv);
  vecs_poly=ones(size(vv));
  for jj=1:npoly,
    vecs_poly(:,jj+1)=vecs_poly(:,jj).*vv;
  end

  vecs=[vecs_poly(:,2:end) vecs_pt];
  assign_bad_timestreams(vecs,tods(j));
  timestreams(j).map=zeros(size(vecs,2),numel(rows));
end

  