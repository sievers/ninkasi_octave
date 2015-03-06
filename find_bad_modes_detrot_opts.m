function[vecs]=find_bad_modes_detrot_opts(tod,opts)
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


is_rot=are_dets_rotated(tod);

rotate_data_detpairs_c (tod);

datft=get_data_fft_c(tod);
n=get_tod_ndata(tod);


dt=get_tod_dt(tod);
nu=(0:(size(datft,1)-1))'/(n*dt);
  
vecs1=find_bad_modes_block(datft(:,~is_rot),nu,freqs,eig_thresh,skip_mean);
vecs2=find_bad_modes_block(datft(:,is_rot),nu,freqs,eig_thresh,false);

if isempty(vecs1)|isempty(vecs2)
  crap=get_tod_name(tod);
  while iscell(crap)
    crap=crap{1};
  end
  disp('we have a problem in find_bad_modes_detrot_opts: ')
  disp(crap)
end
vecs=zeros(length(is_rot),size(vecs1,2)+size(vecs2,2));
vecs(~is_rot,1:size(vecs1,2))=vecs1;
vecs(is_rot,size(vecs1,2)+1:end)=vecs2;
%get raw data back to where it was
rotate_data_detpairs_c (tod);
