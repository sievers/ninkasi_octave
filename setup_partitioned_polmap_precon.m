function[map]=setup_partitioned_polmap_precon(map)
mapptr=allocate_ninkasi_skymap(0.1,0,1,0,1);
p=map.partition;
set_skymap_cea_predef_c(mapptr,p.global_radelt,p.global_decdelt,p.global_rapix,p.global_decpix,p.global_pv,p.my_patch(2)-p.my_patch(1)+1,p.my_patch(4)-p.my_patch(3)+1);

polstate=get_map_polstate_c(map.mapptr);
set_map_polstate(mapptr,polstate);


map.mapptr=mapptr;

tmp=skymap2octave(mapptr);
mm=map.map;
octave2skymap(map.map,map.mapptr);

