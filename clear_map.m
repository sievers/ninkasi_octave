function[map]=clear_map(map)
if isstruct(map)
  clear_map_c(map.mapptr);
else
  clear_map_c(map);
end