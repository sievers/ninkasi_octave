function[filt_org]=debutterworth_octave(tod,debutter,varargin)
if ~exist('debutter')
  debutter=true;
end


dat=get_tod_data(tod);
n=size(dat,1);
dt=get_tod_dt(tod);
filt=butterworth_filter(n,dt,varargin{:});
filt_org=filt;
if debutter
  filt=1./filt;
end


if (1)
  mdisp('doing fast filt');
  apply_complex_filter_to_data(tod,filt);
  return
end


dat=fft_r2c_octave(dat);

filt=filt(1:size(dat,1));
for j=1:size(dat,2),
  dat(:,j)=dat(:,j).*filt;
end
dat=fft_c2r(dat,iseven(n));
push_tod_data(dat,tod);
  