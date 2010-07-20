function[mm]=mapset2mapset_octave(tods,mm,map)
map2=make_map_copy(map);
clear_map(map2);
octave2skymap(mm,map);
for j=1:length(tods),
  mytod=tods(j);
  allocate_tod_storage(mytod);
  map2tod(map,mytod);
  tod2map(mytod,map2);
  free_tod_storage(mytod);
end


mm=skymap2octave(map2);
destroy_map(map2);
