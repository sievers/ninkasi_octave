function[map]=make_map_copy(map)
if isstruct(map)
  map.mapptr=make_map_copy_c(map.mapptr);
else
  map=make_map_copy_c(map);
end