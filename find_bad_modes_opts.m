function[vecs,mat]=find_bad_modes_opts(tod,opts)
mdisp(['welcome to find_bad_modes_opts']);

if ~exist('opts')
  opts.fwee=0;
end

freqs=get_struct_mem(opts,'mode_freqs',[0.25 4 1000]);
eig_thresh=get_struct_mem(opts,'mode_thresh',4 );
skip_mean=get_struct_mem(opts,'mode_skip_mean',false);


nn=length(eig_thresh);
for j=nn+1:length(freqs)-1,
  eig_thresh(j)=eig_thresh(nn);
end


%data=get_tod_data(tod);
%n=size(data,1);
%datft=fft_omp_r2c_many(data);
%clear data;

datft=get_data_fft_c(tod);

n=get_tod_ndata(tod);

dt=get_tod_dt(tod);
nu=(0:(size(datft,1)-1))'/(n*dt);


vecs=find_bad_modes_block(datft,nu,freqs,eig_thresh,skip_mean);

