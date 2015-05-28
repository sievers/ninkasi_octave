function[value]=octave2skymap(map,mapptr)
if isstruct(map)
  if ~isfield(map,'map')
    error('need a map inside the struct in octave2skymap');
  end
  if ~isfield(map,'mapptr')
    error('need a mapptr inside the struct in octave2skymap');
  end
  if isfield(map,'partition')
    my_patch=get_local_map_piece(map);
    octave2skymap_c(my_patch,map.mapptr);
  else
    octave2skymap_c(map.map,map.mapptr);
  end
else
  assert(class(mapptr)=='int64');
  octave2skymap_c(map,mapptr);
end

