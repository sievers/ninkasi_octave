function[bands]=set_tod_noise_bands_projvecs(tod,myopts)
if (does_tod_have_demod(tod))
  myf=get_struct_mem(myopts,'demod_noise_func',@set_tod_noise_bands_projvecs_demod);
  if isa(myf,'function_handle')
    feval(myf,tod,myopts);
  else
    assert(isa(myf,'string'));
    to_exec=[myf '(tod,myopts);'];
    mdisp(to_exec);
    eval(to_exec);
  end
  %set_tod_noise_bands_projvecs_demod(tod,myopts);
  return
end

mdisp('setting projvec noise.');


do_rot=get_struct_mem(myopts,'rotate_detpairs',false);

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


find_modes_new=get_struct_mem(myopts,'find_modes_new',false);


if (find_modes_new)
  %this is the usual path.  Just finding the array patterns here
  if do_rot,
    vecs=find_bad_modes_detrot_opts(tod,myopts);
  else
    vecs=find_bad_modes_opts(tod,myopts);
  end
  to_keep=max(abs(vecs))<det_mode_thresh;  %find modes that are too single-detector
  if min(to_keep)==0,
    %whos vecs
    vecs=vecs(:,to_keep);
    mdisp(['Dropping ' num2str(sum(to_keep==0)) ' modes for being single-detector.']);
    %whos vecs
  end
  nvecs=size(vecs,2);
  if do_rot
    rotate_data_detpairs_c (tod);
  end
  dataft=get_data_fft_c(tod);
  if do_rot
    rotate_data_detpairs_c (tod);
  end
  
  %whos vecs
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

%nband
ncut=0;

%now loop over frequency bins and get the noise properties in each
for j=nband:-1:1,  
  myblock=dataft(ibands(j)+1:ibands(j+1),:);
  %fit modes to the data, measure mode noise.
  %This implementation relies on the modes being orthogonal as they are out of the mode finder
  block_amps=myblock*vecs;
  big_block_amps(j,:)=mean(abs(block_amps).^2)/n;

  if skip_corrnoise_sub
    %mdisp('not subtracting correlated noise from data in set_tod_noise_bands_projvecs.m.');
    myblock_clean=myblock;
  else
    %The usual path - remove the correlated noise from the block
    myblock_clean=myblock-block_amps*vecs';
  end
  %This is the detector noise part of the model
  big_det_noise_amps(j,:)=mean(abs(myblock_clean).^2)/n;
  %if a detector looks suspiciously sensitive, force it to look like the median
  if (det_weight_thresh>0)
    medval=median(big_det_noise_amps(j,:));
    ind=big_det_noise_amps(j,:)>det_weight_thresh*medval;
    big_det_noise_amps(j,ind)=medval;
    ncut=ncut+sum(ind);
  end


  if (0)
    vecs_use=vecs;
    for k=1:nvecs,      
      try
        myfac=sqrt(big_block_amps(j,k));
      catch
        whos
        nvecs
        error(['broken tod is ' get_tod_name(tod)]);
      end
      
      vecs_use(:,k)=vecs(:,k)*myfac*mode_scale_fac;
    end;    
  else
    %the usual path - rather than store eigenmode and amplitudes separately, scale each mode by the sqrt of its
    %amplitude so the correlated noise just becomes VV' instead of V*lambda*V'

    vecs_use=vecs;for k=1:nvecs, vecs_use(:,k)=vecs(:,k)*sqrt(big_block_amps(j,k))*mode_scale_fac;end;
  end
  %Tell the C part of the code what we found
  set_oneband_tod_noise_banded_projvec(tod,j,scale_facs(j)./big_det_noise_amps(j,:),vecs_use);

end
if ncut>0
  mdisp(['Capped ' num2str(100*ncut/numel(big_det_noise_amps)) ' percent of detector weights.']);
end
