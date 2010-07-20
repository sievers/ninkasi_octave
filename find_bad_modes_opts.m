function[vecs]=find_bad_modes_opts(tod,opts)
if ~exist('opts')
  opts.fwee=0;
end

freqs=get_struct_mem(opts,'mode_freqs',[0.25 4 1000]);
eig_thresh=get_struct_mem(opts,'mode_thresh',4 );
nn=length(eig_thresh);
for j=nn+1:length(freqs)-1,
  eig_thresh(j)=eig_thresh(nn);
end

data=get_tod_data(tod);
n=size(data,1);
datft=fft_omp_r2c_many(data);
clear data;
dt=get_tod_dt(tod);
nu=(0:(size(datft,1)-1))'/(n*dt);

vecs=ones(size(datft,2),1);  %explicitly project the mean
for j=1:size(vecs,2),
  vecs(:,j)=vecs(:,j)/sqrt(sum(vecs(:,j).^2));
end

for j=1:length(freqs)-1,
  ind=(nu>freqs(j))&(nu<=freqs(j+1));
  crud=datft(ind,:);
  mat=real(crud'*crud);
  mat=mat+mat';
  mat=project_vecs_from_mat(mat,vecs);
  mat=0.5*(mat+mat');
  [vv,ee]=eig(mat);
  ee=diag(ee);
  if eig_thresh(j)>0,
    val=eig_thresh(j)*median(ee);
    ind=ee>val;
  else
    ind=false(size(ee));
    ind(end+eig_thresh+1:end)=true;
  end
  vecs=[vecs vv(:,ind)];
    
  mdisp(sprintf('cutting %d out of %d modes on band %d\n',sum(ind),length(ind),j));
  clear crud
end

