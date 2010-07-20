function[x]=make_signal_only_maps_fcn()
more off
addpath /cita/d/raid-sievers/sievers/matlab
mpi_init;
myid=mpi_comm_rank+1;
nproc=mpi_comm_size;
format short g
dr='/mnt/raid-cita/sievers/act/hilton_sims/fullTODs_skyMap_NewPointingSignalOnly0Modes/';
tod_names={
'1221390748.1221390761.ar1'
'1221993788.1221993903.ar1'
'1222252186.1222252283.ar1'
'1223630848.1223630964.ar1'
'1223700209.1234326370.ar1'
'1223717251.1234327922.ar1'
'1223975543.1223975689.ar1'
'1224130828.1224130940.ar1'
'1224217241.1224217412.ar1'
'1224475437.1224475573.ar1'
'1224578536.1224578653.ar1'
%'1224923283.1224923419.ar1'
'1225009631.1225009755.ar1'
'1225078521.1225078640.ar1'
'1225267927.1225268037.ar1'
'1225354329.1225354462.ar1'
'1225699022.1234639633.ar1'
'1225768004.1234649611.ar1'
'1226112714.1226112873.ar1'
'1226285106.1226357005.ar1'
'1226302013.1226348217.ar1'
'1226560306.1226619322.ar1'
'1226646707.1226646820.ar1'
'1226991401.1226991516.ar1'
'1227491393.1227491515.ar1'
'1227939086.1227939191.ar1'
'1228025486.1228025589.ar1'
'1228283781.1228283887.ar1'
'1228783817.1228783833.ar1'
};
for j=1:length(tod_names),
  tod_names{j}=[dr tod_names{j}];
end
[tods_org,lims]=read_all_tod_headers(tod_names,1); pixsize=30/60/60*pi/180;

if (nproc==1)
  tods=tods_org;
else
  tods=tods_org(myid:nproc:length(tods_org));
end

ntod=length(tods);
for j=1:ntod,
  purge_cut_detectors(tods(j));
end

map=allocate_ninkasi_skymap(pixsize,lims(1)-1e-3,lims(2)+1e-3, lims(3)-1e-3,lims(4)+1e-3);
mapset.skymap.mapptr=map;
mapset.skymap.map=skymap2octave(mapset.skymap.mapptr);


wt=make_map_copy(map);
make_weightmap_octave(tods,wt);
weight=skymap2octave(wt);
weight=mpi_allreduce(weight);
destroy_map(wt);
precon.skymap.map=weight;


[b,meds]=create_initial_mapset_octave(tods,mapset,'');

x=clear_mapset(mapset,true);
[x,best]=run_pcg_corrnoise_precon_octave(tods,x,b,precon,@apply_precon,'','maxiter',10);
