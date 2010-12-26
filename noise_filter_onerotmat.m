function[value]=noise_filter_onerotmat(tods)
for j=1:length(tods),
  dat_ft=get_data_fft_c(tods(j));
  edges=pull_noise_bin_edges_c(tods(j));
  rotmat=pull_noise_rotmat(tods(j));
  
  dat_ft=dat_ft*rotmat;
  ndata=get_tod_ndata(tods(j));
  dt=get_tod_dt(tods(j));
  dnu=1/(ndata*dt);
  nuvec=(0:(ndata-1))'*dnu;
  
