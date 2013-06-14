function[value]=find_impulse_response(tods,mapset,xpix,ypix,varargin)
if ~isfield(mapset,'skymap')
  error('Need skymap in find_impulse_response.');
end

mapset_in=clear_mapset(mapset,true);
mapset_in.skymap.map(xpix,ypix)=1;
octave2skymap(mapset_in.skymap);

mapset_out=mapset2mapset_corrnoise_octave(tods,mapset_in,varargin{:});

