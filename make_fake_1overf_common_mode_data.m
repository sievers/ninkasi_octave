function[dat,cm]=make_fake_1overf_common_mode_data(tod,varargin)
if numel(varargin)==1
  myopts=varargin{1};
else
  if numel(varargin)==0
    myopts.asdfaadfkgjahfgasdadg='asdfasdfasdfasdf';
  else
    myopts=varargin2opts(varargin{:});
  end
end
if isempty(tod)
  ndata=get_struct_mem(myopts,'sim_ndata',180000);
  ndet=get_struct_mem(myopts,'sim_ndet',512);
  dt=get_struct_mem(myopts,'sim_dt',1/200);
else
  ndata=get_tod_ndata(tod);
  ndet=get_tod_ndet(tod);
  dt=get_tod_dt(tod);
end
create_in_tod=get_struct_mem(myopts,'create_in_tod',~isempty(tod));
det_knee=get_struct_mem(myopts,'det_knee',0.2);
det_slope=get_struct_mem(myopts,'det_slope',-1.0);
common_knee=get_struct_mem(myopts,'common_knee',2.0);
common_slope=get_struct_mem(myopts,'common_slope',-1.5);
white_level=get_struct_mem(myopts,'white_level',400);
tic;
dat=randn(ndata,ndet)*white_level/sqrt(dt);

if (create_in_tod) 
  nn=ceil((ndata+0.2)/2);
  nuvec=(0:(nn-1))';nuvec=nuvec*1/(dt*ndata);
  nuvec(1)=0.5*nuvec(2);
  scale_vec=1+(nuvec/det_knee).^det_slope;

  push_tod_data(dat,tod);
  clear dat;
  apply_real_filter_to_data(tod,scale_vec);
  if (common_knee>0)
    vec=randn(ndata,1)*white_level/sqrt(dt);
    scale_vec=(nuvec/common_knee).^common_slope;
    vecft=fft_r2c(vec);
    vecft=vecft.*scale_vec;
    cm=fft_c2r(vecft);
    add_vector_to_tod_data(tod,cm);
  end
  if (nargout>0)
    dat=get_tod_data(tod);
  end

  return
end



datft=fft(dat);

clear dat;
nn=ceil((ndata+0.2)/2);
nuvec=(0:(nn-1))';nuvec=nuvec*1/(dt*ndata);
nuvec(1)=0.5*nuvec(2);
if iseven(ndata)
  nuvec=[nuvec;flipud(nuvec(2:end-1))];
else
  nuvec=[nuvec;flipud(nuvec(2:end))];
end
assert(length(nuvec)==ndata);

scale_vec=1+(nuvec/det_knee).^det_slope;

for j=1:ndet
  datft(:,j)=datft(:,j).*scale_vec;
end

dat=real(ifft(datft));clear datft;

if (common_knee>0)
  vec=randn(ndata,1)*white_level/sqrt(dt);
  scale_vec=(nuvec/common_knee).^common_slope;
  vecft=fft(vec);
  vecft=vecft.*scale_vec;
  cm=real(ifft(vecft));
  dat=dat+repmat(cm,[1 ndet]);
end



