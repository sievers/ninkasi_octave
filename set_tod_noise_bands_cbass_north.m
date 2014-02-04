function[value]=set_tod_noise_bands_cbass_north(tod,myopts)
mdisp('greetings from set_tod_noise_bands_cbass_north');

bands=myopts.bands;
if ischar(bands),
  if strcmp(bands,'auto')
    bands=get_band_edges_auto(tod,myopts);
  end
end
nband=numel(bands)-1;
scale_facs=get_struct_mem(myopts,'noise_scale_facs',1);
if numel(scale_facs)<nband,
  scale_facs(end+1:nband)=1;
end
dataft=get_data_fft_c(tod);

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
  big_det_noise_amps(j,:)=mean(abs(myblock).^2)/n;
  set_oneband_tod_noise_banded_projvec(tod,j,scale_facs(j)./big_det_noise_amps(j,:),[]);
end

