function[precon,wt_inv,xx,yy]=setup_sigurd_precon(mapset,tods,weight,varargin)
if numel(varargin)>1,
  myopts=varargin2opts(varargin{:});
else
  myopts=varargin{1};
end

if ~isfield(mapset,'skymap')
  warning('no skymap found in mapset in setup_sigurd_precon');
  return;
end

mapset=clear_mapset(mapset);
[xx,yy]=find_weight_map_center(weight,ceil(0.05*min(size(weight))));
[xx yy weight(xx,yy)]

myopts.check_empty=true;  %since most TODs won't hit the delta function, we can save time here.

clear_map(mapset.skymap.mapptr);
mapset.skymap.map(xx,yy)=1;

new_mapset=mapset2mapset_corrnoise_octave(tods,mapset,myopts);
new_map=new_mapset.skymap.map;
%destroy_map(new_mapset.skymap.mapptr); %appear not to need this - map is evidently not getting copied
clear new_mapset;

wt_inv=0*weight;wt_inv(weight>0)=1./weight(weight>0);
new_map=new_map.*wt_inv;
new_map=circshift(new_map,[1-xx 1-yy]);
new_map=new_map/max(max(abs(new_map)));



[crap,mapft]=fast_2d_convolve(new_map,new_map);
%precon.raw_impulse=new_map;
thresh=1e-7;
ii=abs(mapft)<thresh*max(max(abs(mapft)));mdisp(['zeroing ' num2str(sum(sum(ii))) ' fft pixels.'])
mapft=1./mapft;mapft(ii)=0;

precon.impulse=mapft;
clear mapft;
precon.wtroot=sqrt(wt_inv);
if nargout==1,
  clear wt_inv;
end
