function[asdf]=set_tod_noise_bands_projvecs_demod(tod,myopts)
mdisp(['welcome to demod noise'])
det_mode_thresh=get_struct_mem(myopts,'det_mode_thresh',1);
eig_thresh=get_struct_mem(myopts,'mode_thresh',4*4 );
mode_scale_fac=get_struct_mem(myopts,'noise_mode_scale_fac',1.0);  %scale the correlated noise power spectra by this factor
scale_facs=get_struct_mem(myopts,'noise_scale_facs',1);

bands=myopts.bands; %if not here, we will crash.
if ischar(bands),
  if strcmp(bands,'auto'),
    bands=get_band_edges_auto(tod,myopts);
  end
end


dt=get_tod_dt(tod);
n=get_tod_ndata(tod);




demodulate_data_c(tod);
cmat=return_data_from_demod_c(tod);

ibands=ceil(bands*dt*n);
ibands=unique(ibands);  %in case the spacing is too small for some TODs
ibands=ibands(ibands<size(cmat,1)); %now make sure bands go exactly to the end, but not past the demod data
ibands(end+1)=size(cmat,1);
if ibands(1)<1
  ibands(1)=1;
end

if ibands(1)>1,
  ibands=[1 ibands];
end
allocate_tod_noise_banded_projvec(tod,ibands);
%make a variance-normalized dataset, otherwise the eigenvalues will be dominated by T
myvars=mean(abs(cmat).^2);
mysigs=sqrt(myvars);
cmat2=cmat;
for j=1:length(ibands)-1,
  cmat2(:,j)=cmat2(:,j)/mysigs(j);
end
%transform back into a real-space representation
dd=real(ifft([cmat2;flipud(conj(cmat2(2:end,:)))]));
clear cmat2;
mycov=dd'*dd;mycov=0.5*(mycov+mycov');  %make sure the covariance matrix is symmetric
[vv,ee]=eig(mycov);ee=diag(ee);

%find the modes with big eigenvalues
ii=ee>eig_thresh(1)*median(ee);
bad_modes=vv(:,ii);
mdisp(['starting with ' num2str(size(bad_modes,2)) ' bad modes']);
mymax=max(abs(bad_modes));
keep_vec=mymax<det_mode_thresh;
bad_modes=bad_modes(:,keep_vec);
mdisp(['after chopping single-detector modes, down to ' num2str(size(bad_modes,2))]);

%now transform the bad mode array patterns back into covariance space rather than correlation space
ndet=size(cmat,2);
for j=1:ndet,
  bad_modes(j,:)=bad_modes(j,:)*mysigs(j);
end


%now find the correlated noise timestreams and the cleaned detector timestreams
clear dd
dd=real(ifft([cmat;flipud(conj(cmat(2:end,:)))]));
bad_ts=dd*bad_modes;
fitp=inv(bad_ts'*bad_ts)*(bad_ts'*dd);
dd_clean=dd-bad_ts*fitp;
badft=fft(bad_ts);badft=badft(1:ceil(end/2),:);assert(length(badft)==length(cmat));
cleanft=fft(dd_clean);cleanft=cleanft(1:ceil(end/2),:);assert(length(cleanft)==length(cmat));


%now go take RMSs of detectors and correlated noise in bins
nband=length(ibands)-1;

for j=nband:-1:1,
  myblock=cleanft(ibands(j)+1:ibands(j+1),:);
  block_amps=badft(ibands(j)+1:ibands(j+1),:);
  big_block_amps(j,:)=mean(abs(block_amps).^2)/n;
  big_det_noise_amps(j,:)=mean(abs(myblock).^2)/n;
  vecs_use=bad_modes;nvecs=size(bad_modes,2);
  for k=1:nvecs, vecs_use(:,k)=bad_modes(:,k)*sqrt(big_block_amps(j,k))*mode_scale_fac;end;
  set_oneband_tod_noise_banded_projvec(tod,j,scale_facs(j)./big_det_noise_amps(j,:),vecs_use);
end
free_demod_c(tod);