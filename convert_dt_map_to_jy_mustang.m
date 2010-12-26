function[map_jy]=convert_dt_map_to_jy_mustang(map,pixsize,varargin)


unnorm=get_keyval_default('unnorm',true,varargin{:});  %in case you only want to smear
convert_to_jy=get_keyval_default('jy',true,varargin{:});  %false if you want to leave in original units

beam_fwhm=get_keyval_default('fwhm',[9 38],varargin{:});
if numel(beam_fwhm)==1
  amps=1;
else
  amps=[0.97 0.03];
end
assert(numel(amps)==numel(beam_fwhm));
assert(abs(sum(amps)-1)<1e-14);  %make sure beam adds to one


amps=get_keyval_default('amps',amps,varargin{:});
freq=get_keyval_default('freq',90,varargin{:});

jy=jy_from_cmb (freq,1);
map_jy=0;
if pixsize<1e-4
  pixsize=pixsize*180/pi*3600;
  disp(['Converting suspected radians into degrees, pixsize is now ' num2str(pixsize) ' arcseconds.']);
end

if convert_to_jy
  fac=jy*pixsize^2;
else
  fac=1;
end

for j=1:length(amps),
  nsmooth=beam_fwhm(j)/pixsize/sqrt(8*log(2));  
  map_jy=map_jy+smooth_image(map*fac,nsmooth,'unnorm',unnorm)*amps(j);
end


