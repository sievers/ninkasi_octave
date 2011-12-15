function[new_mapset]=copy_mapset(mapset)
new_mapset=mapset;
if isfield(mapset,'skymap')
  new_mapset.skymap.mapptr=make_map_copy(mapset.skymap.mapptr);
end


    
