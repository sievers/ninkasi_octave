function[map2]=mapset2mapset(tods,map)
map2=make_map_copy(map);
clear_map(map2);
for j=1:length(tods),
  mytod=tods(j);
  allocate_tod_storage(mytod);
  map2tod(map,mytod);
  tod2map(mytod,map2);
  free_tod_storage(mytod);

end
