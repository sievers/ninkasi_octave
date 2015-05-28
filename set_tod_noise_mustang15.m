function[bands,fitp_cm,big_fitp,dat_clean,cm2]=set_tod_noise_mustang15(tod,myopts)

mdisp('setting projvec noise.');


do_rot=get_struct_mem(myopts,'rotate_detpairs',false);
mustang_min_freq=get_struct_mem(myopts,'mustang_min_freq',0.01);
mustang_cm_freq=get_struct_mem(myopts,'mustang_cm_freq',5);%don't worry about fitting common mode above this frequency
mustang_cut_freq=get_struct_mem(myopts,'mustang_cut_freq',45); %don't worry about data above this frequency


bands=myopts.bands; %if not here, we will crash.
if ischar(bands),
  if strcmp(bands,'auto'),
    bands=get_band_edges_auto(tod,myopts);
  end
end





nband=numel(bands)-1;
scale_facs=get_struct_mem(myopts,'noise_scale_facs',1); 
det_mode_thresh=get_struct_mem(myopts,'det_mode_thresh',1);  %if a mode has a single amplitude on a detector greater than this, assume 
                                                             %that essentially a single-detector mode has slipped through, so don't include it.
det_weight_thresh=get_struct_mem(myopts,'det_weight_thresh',0);  %if a detector in a band has a weight bigger than this times the median, assume something is screwy, replace by median
                                                                 %ignore if it's <=0
mode_scale_fac=get_struct_mem(myopts,'noise_mode_scale_fac',1.0);  %scale the correlated noise power spectra by this factor
skip_corrnoise_sub=get_struct_mem(myopts,'skip_corrnoise_sub',false);


if mode_scale_fac ~=1,
  mdisp(['going to scale correlated noise spectra by ' num2str(mode_scale_fac)]);
end
mdisp('greetings from projvecs')

if numel(scale_facs)<nband,
  scale_facs(end+1:nband)=1;
end



dat_org=get_tod_data(tod);

tt=get_tod_tvec(tod);
dnu=1./(tt(end)-tt(1));
nus=dnu*(0:length(dat_org)-1)';
filt=ones(size(nus));
filt(nus>mustang_cm_freq)=1e-3;
apply_filter_to_tod(tod,filt);
dat=get_tod_data(tod);


cm1=mean(dat')';
mat=legendre_mat(length(dat),3);mat=[cm1 cm1.*mat(:,2) mat];


fitp=inv(mat'*mat)*(mat'*dat);

pp=fitp;pp(1:2,:)=0;
dd=dat-mat*pp;
for jj=1:size(dd,2),
  dd(:,jj)=dd(:,jj)/fitp(1,jj);
end
cm2=mean(dd')';
mat=legendre_mat(length(dat),3);mat=[cm2 cm2.*mat(:,2) mat];
fitp=inv(mat'*mat)*(mat'*dat);
vecs=fitp(1,:)'; %vecs=vecs/sqrt(sum(vecs.^2));

dat_clean=dat-mat*fitp;
dat_use=dat_org-mat*fitp;


nvecs=size(vecs,2);
push_tod_data(dat_clean,tod);
dataft=get_data_fft_c(tod);

ii=(nus>mustang_min_freq)&(nus<mustang_cm_freq);
global mydat;
cmft=fft(cm2);
mydat=abs(cmft(ii)).^2;
global mynus;
mynus=nus(ii);
fac=median(abs(mydat));mydat=mydat/fac;fitp=fminunc(@mustang_chisq,[1 0 -2]);fitp(1:2)=fitp(1:2)*fac;
fitp_cm=fitp;
big_fitp=zeros(3,size(dat,2));
for j=1:size(dat,2),
  mydat=abs(dataft(ii,j)).^2;
  fac=median(abs(mydat));mydat=mydat/fac;fitp=fminunc(@mustang_chisq,[1 0 -1]);fitp(1:2)=fitp(1:2)*fac;
  big_fitp(:,j)=fitp;
end

%this is the original, unfiltered data with the common mode subtracted
push_tod_data(dat_use,tod);
dataft=get_data_fft_c(tod);

%now restore what we had in there originally...
push_tod_data(dat_org,tod);


dt=get_tod_dt(tod);
n=get_tod_ndata(tod);
nuvec=(0:size(dataft,1)-1)';
nuvec=nuvec/(dt*n);

ibands=round(bands*dt*n);
ibands(end)=size(dataft,1);  %make sure noise bands go exactly to the end of the data
ibands(1)=1;
allocate_tod_noise_banded_projvec(tod,ibands);

%nband
ncut=0;

%now loop over frequency bins and get the noise properties in each
for j=nband:-1:1,  
  if bands(j)<=0
    nu_eff=0.5*bands(j+1);
  else
    nu_eff=0.5*(bands(j)+bands(j+1));
  end

  %the data have already had a common mode subtracted
  myblock_clean=dataft(ibands(j)+1:ibands(j+1),:);

  %fit modes to the data, measure mode noise.
  %This implementation relies on the modes being orthogonal as they are out of the mode finder
  
  %block_amps=myblock*vecs;

  
  %ignore the offset part of the common mode power spectrum as it's probably due to detector noise
  block_amps=fitp_cm(2)*nu_eff^fitp_cm(3);
  %block_amps=fitp_cm(1)+fitp_cm(2)*nu_eff^fitp_cm(3);

  big_block_amps(j,:)=mean(abs(block_amps).^2)/n;

  %This is the detector noise part of the model
  big_det_noise_amps(j,:)=mean(abs(myblock_clean).^2)/n;
  %at low frequencies, rely on the power law fit
  if (nu_eff<mustang_cm_freq)
    tmp=big_fitp(1,:)+big_fitp(2,:).*nu_eff.^big_fitp(3,:);
    big_det_noise_amps(j,:)=(tmp')/n;
  end
  if (nu_eff>mustang_cut_freq)
    big_det_noise_amps(j,:)=big_det_noise_amps(j,:)*1e6;
  end

  %if a detector looks suspiciously sensitive, force it to look like the median
  if (det_weight_thresh>0)
    medval=median(big_det_noise_amps(j,:));
    ind=big_det_noise_amps(j,:)>det_weight_thresh*medval;
    big_det_noise_amps(j,ind)=medval;
    ncut=ncut+sum(ind);
  end

  %the usual path - rather than store eigenmode and amplitudes separately, scale each mode by the sqrt of its
  %amplitude so the correlated noise just becomes VV' instead of V*lambda*V'
  vecs_use=vecs;for k=1:nvecs, vecs_use(:,k)=vecs(:,k)*sqrt(big_block_amps(j,k))*mode_scale_fac;end;



  %Tell the C part of the code what we found
  set_oneband_tod_noise_banded_projvec(tod,j,scale_facs(j)./big_det_noise_amps(j,:),vecs_use);

end
if ncut>0
  mdisp(['Capped ' num2str(100*ncut/numel(big_det_noise_amps)) ' percent of detector weights.']);
end
