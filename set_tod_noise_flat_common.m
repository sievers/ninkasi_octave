function[value]=set_tod_noise_flat_common(tod,thresh_freq)
mdisp('setting noise to be flat with a common mode noise and (optionally) a high-pass value.');

dt=get_tod_dt(tod);
n=get_tod_ndata(tod);
nuvec=(0:size(dataft,1)-1)';
nuvec=nuvec/(dt*n);

if exist('thresh_freq')
  bands=[0 thresh_freq 1e10];
  scale_facs=[0 1];
else
  bands=[0 1e10];
  scale_facs=1;
end


ibands=round(bands*dt*n);

if iseven(n)
  ibands(end)=n/2+1;
else
  ibands(end)=(n+1)/2;
end
ibands(1)=1;
allocate_tod_noise_banded_projvec(tod,ibands);
ndet=get_tod_ndet(tod);
nband=length(bands)-1;
for j=nband:-1:1,
  block_amps=ones(ndet,1);

  myblock=dataft(ibands(j)+1:ibands(j+1),:);
  block_amps=myblock*vecs;
  big_block_amps(j,:)=mean(abs(block_amps).^2)/n;

  myblock_clean=myblock-block_amps*vecs';
  big_det_noise_amps(j,:)=mean(abs(myblock_clean).^2)/n;

  vecs_use=vecs;for k=1:nvecs, vecs_use(:,k)=vecs(:,k)*sqrt(big_block_amps(j,k));end;

  set_oneband_tod_noise_banded_projvec(tod,j,scale_facs(j)./big_det_noise_amps(j,:),vecs_use);

end

