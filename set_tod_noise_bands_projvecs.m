function[value]=set_tod_noise_bands_projvecs(tod,myopts)
mdisp('setting projvec noise.');
bands=myopts.bands; %if not here, we will crash.
nband=numel(bands)-1;
scale_facs=get_struct_mem(myopts,'noise_scale_facs',1); 
if numel(scale_facs)<nband,
  scale_facs(end+1:nband)=1;
end


find_modes_new=get_struct_mem(myopts,'find_modes_new',false);


if (find_modes_new)
  vecs=find_bad_modes_opts(tod,myopts);
  nvecs=size(vecs,2);
  dataft=get_data_fft_c(tod);
else
  nvecs=get_struct_mem(myopts,'nbadmode',10);
  dataft=get_data_fft_c(tod);
  covmat=real(dataft(2:end,:)'*dataft(2:end,:));  %make the real-space covariance matrix, neglecting constant term
  covmat=covmat+covmat';
  [v,e]=eig(covmat);
  vecs=v(:,end-nvecs+1:end);
end

dt=get_tod_dt(tod);
n=get_tod_ndata(tod);
nuvec=(0:size(dataft,1)-1)';
nuvec=nuvec/(dt*n);

ibands=round(bands*dt*n);
ibands(end)=size(dataft,1);  %make sure noise bands go exactly to the end of the data
ibands(1)=1;

allocate_tod_noise_banded_projvec(tod,ibands);


for j=nband:-1:1,
  myblock=dataft(ibands(j)+1:ibands(j+1),:);
  block_amps=myblock*vecs;
  big_block_amps(j,:)=mean(abs(block_amps).^2)/n;

  myblock_clean=myblock-block_amps*vecs';
  big_det_noise_amps(j,:)=mean(abs(myblock_clean).^2)/n;

  vecs_use=vecs;for k=1:nvecs, vecs_use(:,k)=vecs(:,k)*sqrt(big_block_amps(j,k));end;

  set_oneband_tod_noise_banded_projvec(tod,j,scale_facs(j)./big_det_noise_amps(j,:),vecs_use);

end
