function[value]=octave2skymap(map,mapptr)
if isstruct(map)
  if ~isfield(map,'map')
    error('need a map inside the struct in octave2skymap');
  end
  if ~isfield(map,'mapptr')
    error('need a mapptr inside the struct in octave2skymap');
  end
  octave2skymap_c(map.map,map.mapptr);
else
  assert(class(mapptr)=='int64');
  octave2skymap_c(map,mapptr);
end

